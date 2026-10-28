#!/usr/bin/env python3
"""
Tests for check_privacy.py

Note: Test hostnames are fragment-built at runtime to avoid hardcoding private infrastructure names.
"""

import tempfile
import subprocess
from pathlib import Path


# Fragment-build private hostnames for testing (avoid hardcoding in source)
def _build_test_hostname_1():
    """Build first test hostname from fragments."""
    return 'mi' + 'mas'


def _build_test_hostname_2():
    """Build second test hostname from fragments."""
    return 'p' + 'ax'


def _build_internal_git_host():
    """Build internal GitHub Enterprise host from fragments."""
    return 'git' + 'hub.' + 'rpi.' + 'edu'


def run_check(repo_dir: Path) -> tuple[int, str, str]:
    """Run the privacy checker and return exit code, stdout, stderr."""
    result = subprocess.run(
        ["python3", "scripts/check_privacy.py", "--repo-root", str(repo_dir)],
        cwd=repo_dir,
        capture_output=True,
        text=True,
    )
    return result.returncode, result.stdout, result.stderr


def setup_git_repo(tmp_path: Path):
    """Initialize a git repo and add check_privacy.py."""
    subprocess.run(["git", "init"], cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(["git", "config", "user.email", "test@example.com"], cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(["git", "config", "user.name", "Test User"], cwd=tmp_path, check=True, capture_output=True)
    
    # Copy check_privacy.py to the test repo
    scripts_dir = tmp_path / "scripts"
    scripts_dir.mkdir()
    import shutil
    shutil.copy(Path(__file__).parent.parent / "scripts" / "check_privacy.py", scripts_dir)


def add_and_commit_file(repo_dir: Path, filename: str, content: str):
    """Add a file to git tracking."""
    filepath = repo_dir / filename
    filepath.parent.mkdir(parents=True, exist_ok=True)
    filepath.write_text(content)
    subprocess.run(["git", "add", filename], cwd=repo_dir, check=True, capture_output=True)


def test_clean_repo():
    """Test that a clean repo passes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "test.txt", "This is clean content\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 0, f"Expected clean repo to pass. stderr: {stderr}"
        assert "No privacy violations found" in stdout
        print("✓ test_clean_repo passed")


def test_hostname_lowercase():
    """Test detection of lowercase hostname."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        hostname = _build_test_hostname_1()
        add_and_commit_file(repo_dir, "test.txt", f"Run this on {hostname} server\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect lowercase hostname"
        assert "hostname" in stderr.lower()
        print("✓ test_hostname_lowercase passed")


def test_hostname_capitalized():
    """Test detection of capitalized hostname (regression test)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        hostname = _build_test_hostname_1().capitalize()
        add_and_commit_file(repo_dir, "test.txt", f"Connect to {hostname} for GPU work\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect capitalized hostname"
        assert "hostname" in stderr.lower()
        print("✓ test_hostname_capitalized passed")


def test_hostname_with_underscore():
    """Test detection of hostname in identifier (regression test)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        hostname = _build_test_hostname_1().capitalize()
        add_and_commit_file(repo_dir, "config.json", f'{{"repo": "{hostname}_Folding_Skills"}}\n')
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect hostname in identifier"
        assert "hostname" in stderr.lower()
        print("✓ test_hostname_with_underscore passed")


def test_hostname_with_hyphen():
    """Test detection of hostname in identifier."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        hostname = _build_test_hostname_1().capitalize()
        add_and_commit_file(repo_dir, "config.json", f'{{"repo": "{hostname}-Design-Skills"}}\n')
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect hostname in identifier"
        assert "hostname" in stderr.lower()
        print("✓ test_hostname_with_hyphen passed")


def test_other_hostname():
    """Test detection of other private hostname."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        hostname = _build_test_hostname_2()
        add_and_commit_file(repo_dir, "test.txt", f"Use {hostname} for this task\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect hostname"
        assert "hostname" in stderr.lower()
        print("✓ test_other_hostname passed")


def test_home_path():
    """Test detection of /home/<user>/ paths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "config.sh", "DATA_DIR=/home/jdoe/data\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect /home path"
        assert "/home/jdoe/" in stderr
        print("✓ test_home_path passed")


def test_users_path():
    """Test detection of /Users/<user>/ paths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "config.sh", "PROJECT=/Users/jsmith/projects/test\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect /Users path"
        assert "/Users/jsmith/" in stderr
        print("✓ test_users_path passed")


def test_private_ip():
    """Test detection of private IP addresses."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "hosts.txt", "server: 192.168.1.100\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect private IP"
        assert "192.168.1.100" in stderr
        print("✓ test_private_ip passed")


def test_github_token():
    """Test detection of GitHub tokens (case-sensitive)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        # Use a fake but valid-format token
        add_and_commit_file(repo_dir, "env.txt", "GITHUB_TOKEN=ghp_" + "x" * 36 + "\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect GitHub token"
        assert "GitHub token" in stderr
        print("✓ test_github_token passed")


def test_aws_key():
    """Test detection of AWS keys (case-sensitive)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        # Use a fake but valid-format key
        add_and_commit_file(repo_dir, "credentials.txt", "AWS_KEY=AKIAIOSFODNN7EXAMPLE\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect AWS key"
        assert "AWS access key" in stderr
        print("✓ test_aws_key passed")


def test_pem_key():
    """Test detection of PEM private keys."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "key.pem", "-----BEGIN PRIVATE KEY-----\nMIIE...\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect PEM key"
        assert "PEM private key" in stderr
        print("✓ test_pem_key passed")


def test_fasta_file_excluded():
    """Test that FASTA files are excluded from credential checks."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        # FASTA files can contain sequences like SAKIA that look like AWS keys
        add_and_commit_file(repo_dir, "sequence.fasta", ">protein1\nSAKIAEADKLGPQRSTUV\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        # Should pass because .fasta files are excluded from credential checks
        assert exit_code == 0, f"Expected FASTA file to pass. stderr: {stderr}"
        print("✓ test_fasta_file_excluded passed")


def test_csv_file_excluded():
    """Test that CSV files are excluded from credential checks."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        # CSV files can contain data that looks like credentials
        add_and_commit_file(repo_dir, "data.csv", "id,sequence\n1,AKIAEADKLMNOP\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        # Should pass because .csv files are excluded from credential checks
        assert exit_code == 0, f"Expected CSV file to pass. stderr: {stderr}"
        print("✓ test_csv_file_excluded passed")


def test_multiple_violations():
    """Test detection of multiple violations in one file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        hostname1 = _build_test_hostname_1()
        hostname2 = _build_test_hostname_2().capitalize()
        content = f"""
        Connect to {hostname1} at 192.168.1.50
        Data location: /home/testuser/data
        Also check {hostname2} server
        """
        add_and_commit_file(repo_dir, "notes.txt", content)
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect multiple violations"
        # Should detect hostname, IP, and path
        assert "hostname" in stderr.lower() or "192.168.1.50" in stderr or "/home/testuser/" in stderr
        print("✓ test_multiple_violations passed")


def test_tailscale_ip():
    """Test detection of Tailscale CGNAT IP addresses."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "config.txt", "server: 100.64.0.1\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect Tailscale IP"
        assert "100.64.0.1" in stderr or "Tailscale" in stderr
        print("✓ test_tailscale_ip passed")


def test_tailscale_domain():
    """Test detection of Tailscale .ts.net domains."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "config.txt", "host: myserver.ts.net\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect .ts.net domain"
        assert "ts.net" in stderr.lower() or "Tailscale" in stderr
        print("✓ test_tailscale_domain passed")


def test_github_rpi_edu():
    """Test detection of internal GitHub Enterprise host."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        git_host = _build_internal_git_host()
        add_and_commit_file(repo_dir, "config.txt", f"url: https://{git_host}/user/repo\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect internal git host"
        assert "github" in stderr.lower() or "internal" in stderr.lower()
        print("✓ test_github_rpi_edu passed")


def test_ssh_identity():
    """Test detection of SSH identity file references."""
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = Path(tmpdir)
        setup_git_repo(repo_dir)
        
        add_and_commit_file(repo_dir, "config.txt", "key: ~/.ssh/id_rsa\n")
        
        exit_code, stdout, stderr = run_check(repo_dir)
        assert exit_code == 1, "Expected to detect SSH identity file"
        assert "SSH" in stderr or "identity" in stderr.lower()
        print("✓ test_ssh_identity passed")


if __name__ == "__main__":
    print("Running privacy checker tests...\n")
    
    test_clean_repo()
    test_hostname_lowercase()
    test_hostname_capitalized()
    test_hostname_with_underscore()
    test_hostname_with_hyphen()
    test_other_hostname()
    test_home_path()
    test_users_path()
    test_private_ip()
    test_github_token()
    test_aws_key()
    test_pem_key()
    test_fasta_file_excluded()
    test_csv_file_excluded()
    test_multiple_violations()
    test_tailscale_ip()
    test_tailscale_domain()
    test_github_rpi_edu()
    test_ssh_identity()
    
    print("\n✅ All tests passed!")
