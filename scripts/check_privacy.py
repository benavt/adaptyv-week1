#!/usr/bin/env python3
"""
Privacy checker for public repositories.

Scans all tracked files (text, binary, PDF, Office) for private infrastructure details:
- Private hostnames - case-insensitive, including embedded forms
- Absolute user paths (/home/<user>/, /Users/<user>/)
- Private IP addresses (RFC 1918, Tailscale 100.64/10)
- Tailscale domains (.ts.net)
- Internal git hosts
- SSH identity file references
- Common credential patterns (GitHub tokens, AWS keys, PEM keys) - case-sensitive

Scanning methods:
- Text files: Direct content scan
- Binary files: Extract printable ASCII strings (≥4 chars) and scan
- PDF files: Text extraction + embedded attachments
- Office files: Unzip and scan XML content

Exits with code 1 if any violations are found.
"""

import argparse
import ipaddress
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import List, Tuple


# Private hostnames to detect (case-insensitive with word boundaries or underscores/hyphens)
# Fragment-built to avoid hardcoding sensitive infrastructure names
PRIVATE_HOSTNAMES = ['mi' + 'mas', 'p' + 'ax']

# File extensions to exclude from credential pattern checks (legitimate false positives)
# Protein sequence files and CSV data files can contain substrings like SAKIA or AKIAEADK
# that resemble AWS keys but are actually scientific data
CREDENTIAL_CHECK_EXCLUSIONS = {".fasta", ".a3m", ".pdb", ".cif", ".csv"}

# Files to completely skip (contain test data and documentation of violations)
SKIP_FILES = {
    "scripts/check_privacy.py",  # This script contains examples of violations
    "tests/test_check_privacy.py",  # Tests contain example violations
}


def get_tracked_files(repo_root: Path) -> List[Path]:
    """Get all tracked files from git."""
    try:
        result = subprocess.run(
            ["git", "ls-files"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        )
        return [repo_root / line for line in result.stdout.strip().split("\n") if line]
    except subprocess.CalledProcessError as e:
        print(f"Error getting tracked files: {e}", file=sys.stderr)
        sys.exit(1)


def is_text_file(path: Path) -> bool:
    """Check if a file is likely a text file."""
    try:
        with open(path, "rb") as f:
            # Read first 8KB
            chunk = f.read(8192)
            if not chunk:
                return True  # Empty file
            # Check for null bytes (binary indicator)
            return b"\x00" not in chunk
    except Exception:
        return False


def extract_printable_strings(filepath: Path, min_length: int = 4) -> str:
    """
    Extract printable ASCII strings from a binary file (like the 'strings' command).
    Returns extracted strings joined by newlines for pattern matching.
    """
    try:
        with open(filepath, "rb") as f:
            data = f.read()
        
        # Extract sequences of printable ASCII characters (0x20-0x7E)
        strings = []
        current = []
        
        for byte in data:
            if 0x20 <= byte <= 0x7E:  # Printable ASCII
                current.append(chr(byte))
            else:
                if len(current) >= min_length:
                    strings.append(''.join(current))
                current = []
        
        # Don't forget the last string
        if len(current) >= min_length:
            strings.append(''.join(current))
        
        return '\n'.join(strings)
    except Exception:
        return ""


def extract_png_text_chunks(filepath: Path) -> str:
    """
    Parse PNG file and extract text from text-bearing chunks only.
    
    Scans:
    - tEXt: Latin-1 text
    - iTXt: UTF-8 text (optionally compressed)
    - zTXt: Compressed Latin-1 text
    - eXIf: EXIF data (may contain text)
    - Unknown ancillary chunks (lowercase first letter)
    
    Skips:
    - IDAT: Compressed image data (would cause false positives)
    - Other critical chunks (IHDR, PLTE, IEND, etc.)
    """
    import struct
    import zlib
    
    PNG_SIGNATURE = b'\x89PNG\r\n\x1a\n'
    
    # Critical chunks to skip (uppercase first letter = critical)
    SKIP_CHUNKS = {b'IHDR', b'PLTE', b'IDAT', b'IEND', b'cHRM', b'gAMA', 
                   b'iCCP', b'sBIT', b'sRGB', b'bKGD', b'hIST', b'tRNS',
                   b'pHYs', b'sPLT', b'tIME', b'acTL', b'fcTL', b'fdAT'}
    
    try:
        with open(filepath, "rb") as f:
            # Check PNG signature
            signature = f.read(8)
            if signature != PNG_SIGNATURE:
                return ""  # Not a PNG
            
            text_content = []
            
            while True:
                # Read chunk length
                length_bytes = f.read(4)
                if len(length_bytes) < 4:
                    break
                
                length = struct.unpack('>I', length_bytes)[0]
                chunk_type = f.read(4)
                if len(chunk_type) < 4:
                    break
                
                chunk_data = f.read(length)
                crc = f.read(4)  # Read but don't verify CRC
                
                if len(chunk_data) < length:
                    break
                
                # Skip IDAT and other known non-text chunks
                if chunk_type in SKIP_CHUNKS:
                    continue
                
                # Process text chunks
                if chunk_type == b'tEXt':
                    # tEXt: keyword\0text (Latin-1)
                    try:
                        keyword, _, text = chunk_data.partition(b'\x00')
                        text_content.append(text.decode('latin-1', errors='ignore'))
                    except:
                        pass
                
                elif chunk_type == b'iTXt':
                    # iTXt: keyword\0compression_flag\0compression_method\0language_tag\0translated_keyword\0text
                    try:
                        parts = chunk_data.split(b'\x00', 4)
                        if len(parts) >= 5:
                            compression_flag = parts[1][0] if parts[1] else 0
                            text = parts[4]
                            if compression_flag:
                                text = zlib.decompress(text)
                            text_content.append(text.decode('utf-8', errors='ignore'))
                    except:
                        pass
                
                elif chunk_type == b'zTXt':
                    # zTXt: keyword\0compression_method\0compressed_text
                    try:
                        keyword, rest = chunk_data.split(b'\x00', 1)
                        if rest:
                            compression_method = rest[0]
                            compressed_text = rest[1:]
                            text = zlib.decompress(compressed_text)
                            text_content.append(text.decode('latin-1', errors='ignore'))
                    except:
                        pass
                
                elif chunk_type == b'eXIf':
                    # eXIf: EXIF data (may contain text fields)
                    try:
                        text_content.append(chunk_data.decode('latin-1', errors='ignore'))
                    except:
                        pass
                
                else:
                    # Unknown ancillary chunk (lowercase first letter)
                    if chunk_type[0:1].islower():
                        try:
                            text_content.append(chunk_data.decode('latin-1', errors='ignore'))
                        except:
                            pass
            
            return '\n'.join(text_content)
    
    except Exception:
        return ""


def check_hostname_patterns(content: str, line_num: int, filepath: Path) -> List[Tuple[int, str]]:
    """Check for private hostname patterns (case-insensitive)."""
    violations = []
    content_lower = content.lower()
    
    for hostname in PRIVATE_HOSTNAMES:
        # Pattern matches hostname with word boundaries or as part of identifiers
        pattern = r"(?:^|[^a-zA-Z0-9])(" + re.escape(hostname) + r")(?:[^a-zA-Z0-9]|$)"
        if re.search(pattern, content_lower):
            violations.append((line_num, "Private hostname found"))
    
    return violations


def check_user_paths(content: str, line_num: int) -> List[Tuple[int, str]]:
    """Check for absolute user paths."""
    violations = []
    
    # Match /home/<username>/ or /Users/<username>/
    home_pattern = r"/home/([a-zA-Z0-9_-]+)/"
    users_pattern = r"/Users/([a-zA-Z0-9_-]+)/"
    
    for match in re.finditer(home_pattern, content):
        violations.append((line_num, f"Private path found: {match.group(0)}"))
    
    for match in re.finditer(users_pattern, content):
        violations.append((line_num, f"Private path found: {match.group(0)}"))
    
    return violations


def check_private_ips(content: str, line_num: int) -> List[Tuple[int, str]]:
    """Check for private IP addresses (RFC 1918 and Tailscale 100.64/10)."""
    violations = []
    
    # Pattern for IPv4 addresses
    ip_pattern = r"\b(?:[0-9]{1,3}\.){3}[0-9]{1,3}\b"
    
    for match in re.finditer(ip_pattern, content):
        ip_str = match.group(0)
        try:
            ip = ipaddress.ip_address(ip_str)
            # Check RFC 1918 private ranges
            if ip.is_private and not ip.is_loopback:
                violations.append((line_num, f"Private IP address found: {ip_str}"))
            # Check Tailscale CGNAT range (100.64.0.0/10)
            elif ipaddress.ip_address(ip_str) in ipaddress.ip_network("100.64.0.0/10"):
                violations.append((line_num, f"Tailscale IP address found: {ip_str}"))
        except ValueError:
            # Not a valid IP address
            pass
    
    return violations


def check_private_domains(content: str, line_num: int) -> List[Tuple[int, str]]:
    """Check for private domain patterns."""
    violations = []
    
    # Tailscale domains (.ts.net)
    if re.search(r"\b[a-zA-Z0-9\-]+\.ts\.net\b", content, re.IGNORECASE):
        violations.append((line_num, "Tailscale domain (.ts.net) found"))
    
    # Internal GitHub Enterprise host (fragment-built)
    internal_git_host = 'git' + 'hub.' + 'rpi.' + 'edu'
    if re.search(r"\b" + re.escape(internal_git_host) + r"\b", content, re.IGNORECASE):
        violations.append((line_num, "Internal GitHub Enterprise host found"))
    
    return violations


def check_ssh_identities(content: str, line_num: int) -> List[Tuple[int, str]]:
    """Check for SSH identity file references."""
    violations = []
    
    # Common SSH identity file patterns
    ssh_patterns = [
        r"\.ssh/id_rsa",
        r"\.ssh/id_dsa",
        r"\.ssh/id_ecdsa",
        r"\.ssh/id_ed25519",
        r"IdentityFile\s+",
    ]
    
    for pattern in ssh_patterns:
        if re.search(pattern, content, re.IGNORECASE):
            violations.append((line_num, "SSH identity file reference found"))
            break
    
    return violations


def check_credentials(content: str, line_num: int, filepath: Path) -> List[Tuple[int, str]]:
    """Check for credential patterns (case-sensitive).
    
    Excludes protein sequence files and CSV files where false positives are common.
    """
    # Skip credential checks for excluded file types
    if filepath.suffix.lower() in CREDENTIAL_CHECK_EXCLUSIONS:
        return []
    
    violations = []
    
    # PEM private key headers (case-sensitive)
    if "-----BEGIN PRIVATE KEY-----" in content or "-----BEGIN RSA PRIVATE KEY-----" in content:
        violations.append((line_num, "PEM private key header found"))
    
    # GitHub tokens (case-sensitive)
    github_token_pattern = r"\b(ghp_[a-zA-Z0-9]{36}|gho_[a-zA-Z0-9]{36}|ghu_[a-zA-Z0-9]{36}|ghs_[a-zA-Z0-9]{36}|ghr_[a-zA-Z0-9]{36}|github_pat_[a-zA-Z0-9_]{82})\b"
    for match in re.finditer(github_token_pattern, content):
        violations.append((line_num, f"GitHub token found: {match.group(0)[:20]}..."))
    
    # AWS access keys (case-sensitive)
    aws_key_pattern = r"\b(AKIA[0-9A-Z]{16})\b"
    for match in re.finditer(aws_key_pattern, content):
        violations.append((line_num, f"AWS access key found: {match.group(0)[:10]}..."))
    
    return violations


def extract_pdf_text(filepath: Path) -> str:
    """Extract text from PDF using pdftotext if available."""
    if not shutil.which("pdftotext"):
        return ""
    
    try:
        result = subprocess.run(
            ["pdftotext", str(filepath), "-"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.stdout if result.returncode == 0 else ""
    except Exception:
        return ""


def extract_pdf_attachments(filepath: Path, temp_dir: Path) -> List[Path]:
    """Extract embedded PDF attachments using pdfdetach if available."""
    if not shutil.which("pdfdetach"):
        return []
    
    try:
        # First check if there are attachments
        result = subprocess.run(
            ["pdfdetach", "-list", str(filepath)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return []
        
        # Save attachments to temp directory
        subprocess.run(
            ["pdfdetach", "-saveall", "-o", str(temp_dir), str(filepath)],
            capture_output=True,
            timeout=30,
        )
        
        return list(temp_dir.glob("*"))
    except Exception:
        return []


def extract_office_xml(filepath: Path, temp_dir: Path) -> List[Path]:
    """Extract XML files from Office documents (docx, xlsx, pptx)."""
    try:
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            # Extract only XML files
            xml_files = []
            for name in zip_ref.namelist():
                if name.endswith('.xml'):
                    extract_path = temp_dir / name
                    extract_path.parent.mkdir(parents=True, exist_ok=True)
                    extract_path.write_bytes(zip_ref.read(name))
                    xml_files.append(extract_path)
            return xml_files
    except Exception:
        return []


def scan_text_content(content: str, filepath: Path, line_offset: int = 0) -> List[Tuple[int, str]]:
    """Scan text content for all violation patterns."""
    violations = []
    
    for line_num, line in enumerate(content.split('\n'), start=line_offset + 1):
        violations.extend(check_hostname_patterns(line, line_num, filepath))
        violations.extend(check_user_paths(line, line_num))
        violations.extend(check_private_ips(line, line_num))
        violations.extend(check_private_domains(line, line_num))
        violations.extend(check_ssh_identities(line, line_num))
        violations.extend(check_credentials(line, line_num, filepath))
    
    return violations


def scan_file(filepath: Path, repo_root: Path) -> List[Tuple[int, str]]:
    """Scan a file for privacy violations."""
    # Skip files that contain test data or documentation of violations
    rel_path = filepath.relative_to(repo_root)
    if str(rel_path) in SKIP_FILES:
        return []
    
    violations = []
    suffix = filepath.suffix.lower()
    
    # Handle PDF files with text and attachment extraction
    if suffix == '.pdf':
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Extract and scan PDF text
            pdf_text = extract_pdf_text(filepath)
            if pdf_text:
                violations.extend(scan_text_content(pdf_text, filepath))
            
            # Extract and scan PDF attachments
            attachments = extract_pdf_attachments(filepath, temp_path)
            for attachment in attachments:
                if is_text_file(attachment):
                    try:
                        content = attachment.read_text(encoding="utf-8", errors="ignore")
                        violations.extend(scan_text_content(content, filepath))
                    except Exception:
                        pass
        
        return violations
    
    # Handle Office files (docx, xlsx, pptx)
    if suffix in {'.docx', '.xlsx', '.pptx'}:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            xml_files = extract_office_xml(filepath, temp_path)
            
            for xml_file in xml_files:
                try:
                    content = xml_file.read_text(encoding="utf-8", errors="ignore")
                    violations.extend(scan_text_content(content, filepath))
                except Exception:
                    pass
        
        return violations
    
    # Check if file is text or binary
    if is_text_file(filepath):
        # Handle regular text files
        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
                violations = scan_text_content(content, filepath)
        except Exception as e:
            print(f"Warning: Could not read {filepath}: {e}", file=sys.stderr)
    else:
        # Check for PNG by magic bytes
        try:
            with open(filepath, "rb") as f:
                magic = f.read(8)
            
            if magic == b'\x89PNG\r\n\x1a\n':
                # PNG file - parse chunks and scan only text-bearing chunks
                png_text = extract_png_text_chunks(filepath)
                if png_text:
                    violations = scan_text_content(png_text, filepath)
            else:
                # Other binary files - extract printable strings
                strings_content = extract_printable_strings(filepath)
                if strings_content:
                    violations = scan_text_content(strings_content, filepath)
        except Exception as e:
            print(f"Warning: Could not process binary file {filepath}: {e}", file=sys.stderr)
    
    return violations


def main():
    parser = argparse.ArgumentParser(description="Check for private infrastructure details in tracked files")
    parser.add_argument("--repo-root", type=Path, default=Path.cwd(), help="Repository root directory")
    args = parser.parse_args()
    
    repo_root = args.repo_root.resolve()
    if not (repo_root / ".git").exists():
        print(f"Error: {repo_root} is not a git repository", file=sys.stderr)
        sys.exit(1)
    
    print(f"Scanning tracked files in {repo_root}...")
    tracked_files = get_tracked_files(repo_root)
    print(f"Found {len(tracked_files)} tracked files")
    
    all_violations = {}
    for filepath in tracked_files:
        violations = scan_file(filepath, repo_root)
        if violations:
            rel_path = filepath.relative_to(repo_root)
            all_violations[rel_path] = violations
    
    if all_violations:
        print("\n❌ Privacy violations found:\n", file=sys.stderr)
        for filepath, violations in sorted(all_violations.items()):
            print(f"{filepath}:", file=sys.stderr)
            for line_num, message in violations:
                print(f"  Line {line_num}: {message}", file=sys.stderr)
        
        print(f"\nTotal: {sum(len(v) for v in all_violations.values())} violations in {len(all_violations)} files", file=sys.stderr)
        sys.exit(1)
    else:
        print("✓ No privacy violations found")
        sys.exit(0)


if __name__ == "__main__":
    main()
