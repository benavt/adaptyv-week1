#!/usr/bin/env python3
"""
Privacy checker for public repositories.

Scans tracked text files for private infrastructure details that should not appear in public repos:
- Private hostnames (mimas, pax) - case-insensitive, including embedded forms
- Absolute user paths (/home/<user>/, /Users/<user>/)
- Private IP addresses (RFC 1918)
- Common credential patterns (GitHub tokens, AWS keys, PEM keys) - case-sensitive

Exits with code 1 if any violations are found.
"""

import argparse
import ipaddress
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


# Private hostnames to detect (case-insensitive with word boundaries or underscores/hyphens)
PRIVATE_HOSTNAMES = ["mimas", "pax"]

# File extensions to exclude from credential pattern checks (legitimate false positives)
# Protein sequence files and CSV data files can contain substrings like SAKIA or AKIAEADK
# that resemble AWS keys but are actually scientific data
CREDENTIAL_CHECK_EXCLUSIONS = {".fasta", ".a3m", ".pdb", ".cif", ".csv"}


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


def check_hostname_patterns(content: str, line_num: int, filepath: Path) -> List[Tuple[int, str]]:
    """Check for private hostname patterns (case-insensitive)."""
    violations = []
    content_lower = content.lower()
    
    for hostname in PRIVATE_HOSTNAMES:
        # Pattern matches hostname with word boundaries or as part of identifiers
        # e.g., "mimas", "Mimas", "MIMAS", "Mimas_Folding_Skills", "Mimas-Design-Skills"
        pattern = r"(?:^|[^a-zA-Z0-9])(" + re.escape(hostname) + r")(?:[^a-zA-Z0-9]|$)"
        if re.search(pattern, content_lower):
            violations.append((line_num, f"Private hostname '{hostname}' found"))
    
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
    """Check for private IP addresses (RFC 1918)."""
    violations = []
    
    # Pattern for IPv4 addresses
    ip_pattern = r"\b(?:[0-9]{1,3}\.){3}[0-9]{1,3}\b"
    
    for match in re.finditer(ip_pattern, content):
        ip_str = match.group(0)
        try:
            ip = ipaddress.ip_address(ip_str)
            if ip.is_private and not ip.is_loopback:
                violations.append((line_num, f"Private IP address found: {ip_str}"))
        except ValueError:
            # Not a valid IP address
            pass
    
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


def scan_file(filepath: Path, repo_root: Path) -> List[Tuple[int, str]]:
    """Scan a file for privacy violations."""
    if not is_text_file(filepath):
        return []
    
    violations = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line_num, line in enumerate(f, start=1):
                violations.extend(check_hostname_patterns(line, line_num, filepath))
                violations.extend(check_user_paths(line, line_num))
                violations.extend(check_private_ips(line, line_num))
                violations.extend(check_credentials(line, line_num, filepath))
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}", file=sys.stderr)
    
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
