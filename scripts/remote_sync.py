"""Local transport wrapper. Site-specific values come only from environment variables."""

import argparse
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import subprocess
import sys
import tempfile


METADATA = (
    "manifest.json", "plan.json", "requests.jsonl", "requests.csv",
    "jobs.jsonl", "jobs.csv", "results.jsonl", "results.csv",
    "scores.jsonl", "scores.csv", "paired_results.jsonl", "paired_results.csv",
    "summary.json", "collection_manifest.json", "task.json", "state.json",
    "runtime.json",
)
ARTIFACTS = ("structure_path", "confidence_path", "confidence_detail_path", "pae_path")
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15"]


def required(name):
    value = os.environ.get(name, "")
    if not value:
        raise ValueError(f"Set {name} in the local configuration.")
    return value


def remote_path(name):
    value = required(name)
    path = PurePosixPath(value)
    if not path.is_absolute() or path == PurePosixPath("/"):
        raise ValueError(f"{name} must be an absolute non-root directory.")
    if ".." in path.parts or any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError(f"{name} contains an unsupported path component.")
    return path


def local_root(name, default):
    path = Path(os.environ.get(name) or default).expanduser()
    if not path.is_absolute():
        path = Path(os.environ["PROJECT_ROOT"]) / path
    return path.resolve()


def identifier(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise argparse.ArgumentTypeError("Use a run ID starting with a letter/number, then letters/numbers/_.-.")
    return value


def selected_files(index, remote_root, rank=1, engine=None, species=None):
    """Keep complete predictions and sidecars, preserving paths beneath the run root."""
    selected = set()
    with index.open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if not isinstance(row, dict):
                raise ValueError(f"Result index line {number} must be an object.")
            if row.get("status") != "complete" or row.get("rank") != rank:
                continue
            if engine and row.get("engine") != engine:
                continue
            if species and row.get("target_species") != species:
                continue
            if not row.get("structure_path"):
                raise ValueError(f"Complete result on line {number} has no structure_path.")
            for key in ARTIFACTS:
                value = row.get(key)
                if not value:
                    continue
                if not isinstance(value, str) or any(ord(c) < 32 or ord(c) == 127 for c in value):
                    raise ValueError(f"Invalid {key} on result index line {number}.")
                path = PurePosixPath(value)
                if ".." in path.parts:
                    raise ValueError(f"Unsafe {key} on result index line {number}.")
                try:
                    relative = path.relative_to(remote_root)
                except ValueError:
                    raise ValueError(f"{key} on line {number} is outside the configured remote root; use a separate mapping.") from None
                if relative == PurePosixPath("."):
                    raise ValueError(f"{key} on line {number} names the run root, not an artifact.")
                selected.add(relative.as_posix())
    return sorted(selected)


def read_file_list(path):
    """Accept newline or NUL lists, restricted to paths beneath the source root."""
    content = path.read_bytes()
    names = content.split(b"\0") if b"\0" in content else content.splitlines()
    selected = set()
    for raw in names:
        if not raw:
            continue
        value = raw.decode("utf-8")
        relative = PurePosixPath(value)
        if relative.is_absolute() or ".." in relative.parts or relative == PurePosixPath("."):
            raise ValueError("File lists must contain relative artifact paths beneath the source root.")
        if any(ord(c) < 32 or ord(c) == 127 for c in value):
            raise ValueError("File list contains an unsupported control character.")
        selected.add(relative.as_posix())
    return sorted(selected)


class Transport:
    def __init__(self):
        self.host = required("REMOTE")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", self.host):
            raise ValueError("REMOTE must be a short SSH config alias, without username or address.")
        self.protected_args = None

    def ssh(self, command, capture=False):
        return subprocess.run([*SSH, self.host, command], check=True, text=True,
                              stdout=subprocess.PIPE if capture else None)

    def prepare_upload(self, remote):
        quoted = shlex.quote(str(remote))
        # Inputs belong in a writable staging root, never a linked retired run.
        self.ssh(
            f"p={quoted}; while [ \"$p\" != / ]; do "
            "test ! -L \"$p\" || exit 1; p=${p%/*}; "
            "[ -n \"$p\" ] || p=/; done; "
            f"mkdir -p -- {quoted} && test -w {quoted}"
        )

    def operand(self, path):
        if self.protected_args is None:
            help_result = subprocess.run(["rsync", "--help"], capture_output=True, text=True, check=False)
            self.protected_args = "--protect-args" in help_result.stdout + help_result.stderr
        # Modern rsync transmits raw paths with -s. Legacy rsync/openrsync
        # forwards this operand through the remote shell, so quote it there.
        value = str(path)
        return f"{self.host}:{value if self.protected_args else shlex.quote(value)}"

    def sync(self, remote, local, args, push=False, filters=(), files=None):
        operand = self.operand(str(remote).rstrip("/") + "/")
        command = ["rsync", "-a", "--partial", "--itemize-changes", "--stats",
                   "-e", shlex.join(SSH)]
        if self.protected_args:
            command.append("--protect-args")
        if not push:
            command.append("-L")
        if args.dry_run:
            command.append("--dry-run")
        if args.checksum:
            command.append("--checksum")
        command.extend(filters)
        if files is not None:
            command.extend(["-r", "-0", f"--files-from={files}"])
        local_operand = str(local).rstrip("/") + "/"
        command.extend([local_operand, operand] if push else [operand, local_operand])
        subprocess.run(command, check=True)


def destinations(target, run_id):
    local = local_root("LOCAL_RESULT_ROOT", "results")
    if target == "folding":
        return remote_path("REMOTE_FOLDING_ROOT") / run_id / "folding", local / target / run_id
    if target == "campaign":
        return remote_path("REMOTE_CAMPAIGN_ROOT") / run_id, local / "campaigns" / run_id
    if target == "bundle":
        return remote_path("REMOTE_BUNDLE"), local / "bundles" / run_id
    foundry = remote_path("REMOTE_FOUNDRY_ROOT") / run_id
    if target == "rfd3":
        if os.environ.get("REMOTE_RFD3_ROOT"):
            remote = remote_path("REMOTE_RFD3_ROOT") / run_id / "outputs"
        else:
            remote = foundry / "rfd3" / "outputs"
        return remote, local / "foundry" / run_id / "rfd3_compressed"
    return foundry, local / "foundry" / run_id


def parser():
    cli = argparse.ArgumentParser(description="Upload staged inputs or download run artifacts using private local configuration.")
    commands = cli.add_subparsers(dest="command", required=True)
    transfer = argparse.ArgumentParser(add_help=False)
    transfer.add_argument("--dry-run", action="store_true", help="Preview transfers; never mkdir on the remote host.")
    transfer.add_argument("--checksum", action="store_true", help="Compare file contents instead of size/mtime.")
    prepare = commands.add_parser("prepare", help="Create only the writable input staging destination before a first preview.")
    prepare.add_argument("run_id", type=identifier)
    push = commands.add_parser("push", parents=[transfer], help="Upload one staged experiment; no job submission.")
    push.add_argument("run_id", type=identifier)
    pull = commands.add_parser("pull", parents=[transfer], help="Pull a run, campaign, export, or resolved archive bundle.")
    pull.add_argument("target", choices=["folding", "foundry", "campaign", "rfd3", "bundle"])
    pull.add_argument("run_id", type=identifier)
    modes = pull.add_mutually_exclusive_group()
    modes.add_argument("--metadata", action="store_true", help="Download indexes and execution metadata only.")
    modes.add_argument("--selected", action="store_true", help="Use an already downloaded results.jsonl to select models.")
    modes.add_argument("--files-from", type=Path, help="Pull only the relative paths in a private newline or NUL file list.")
    pull.add_argument("--rank", type=int, default=1)
    pull.add_argument("--engine")
    pull.add_argument("--species", help="Exact target_species value from the index.")
    log = commands.add_parser("log", parents=[transfer], help="Pull an ordinary folding Slurm log.")
    log.add_argument("job_id", type=lambda value: value if value.isdecimal() else identifier(value))
    resolve = commands.add_parser("resolve", help="Print the resolved Foundry run path for archive discovery.")
    resolve.add_argument("run_id", type=identifier)
    commands.add_parser("doctor", help="Check noninteractive SSH, remote rsync, and configured read roots.")
    commands.add_parser("run", help="Explicitly execute REMOTE_CMD from the trusted local config.")
    return cli


def main():
    cli = parser()
    args = cli.parse_args()
    transport = Transport()
    if args.command == "prepare":
        remote = remote_path("REMOTE_INPUT_ROOT") / args.run_id
        transport.prepare_upload(remote)
        print("Input staging destination ready.")
    elif args.command == "push":
        local = local_root("LOCAL_INPUT_ROOT", "inputs") / args.run_id
        remote = remote_path("REMOTE_INPUT_ROOT") / args.run_id
        if not local.is_dir() or local.is_symlink():
            raise ValueError("Create a real staged input directory for this run ID first.")
        if not any(local.iterdir()):
            raise ValueError("Staged input directory is empty.")
        if any(path.is_symlink() for path in local.rglob("*")):
            raise ValueError("Materialize local input symlinks before uploading.")
        if not args.dry_run:
            transport.prepare_upload(remote)
        transport.sync(remote, local, args, push=True)
    elif args.command == "pull":
        if (args.metadata or args.selected) and args.target not in ("folding", "campaign"):
            cli.error("--metadata/--selected apply only to folding or campaign.")
        if not args.selected and (args.rank != 1 or args.engine or args.species):
            cli.error("--rank/--engine/--species require --selected.")
        if args.rank < 1:
            cli.error("--rank must be positive.")
        remote, local = destinations(args.target, args.run_id)
        filters = []
        if args.metadata:
            filters = ["--include=*/", *(f"--include={name}" for name in METADATA),
                       "--exclude=*", "--prune-empty-dirs"]
        if args.target == "rfd3":
            filters = ["--include=*/", "--include=*.cif.gz", "--include=*.json",
                       "--exclude=*", "--prune-empty-dirs"]
        if args.selected or args.files_from:
            if args.selected:
                index = local / ("folding/results.jsonl" if args.target == "campaign" else "results.jsonl")
                files = selected_files(index, remote, args.rank, args.engine, args.species)
            else:
                files = read_file_list(args.files_from)
            if not files:
                print("No artifacts selected; nothing transferred.")
                return
            # Private absolute paths never enter a committed selection manifest.
            with tempfile.NamedTemporaryFile(mode="wb") as manifest:
                manifest.write(("\0".join(files) + "\0").encode("utf-8"))
                manifest.flush()
                local.mkdir(parents=True, exist_ok=True)
                transport.sync(remote, local, args, files=manifest.name)
        else:
            local.mkdir(parents=True, exist_ok=True)
            transport.sync(remote, local, args, filters=filters)
    elif args.command == "log":
        local = local_root("LOCAL_LOG_ROOT", "logs")
        local.mkdir(parents=True, exist_ok=True)
        remote = remote_path("REMOTE_LOG_ROOT") / f"folding-{args.job_id}.out"
        # Logs are files, so this transfer intentionally has no source slash.
        operand = transport.operand(remote)
        command = ["rsync", "-aL", "--partial", "--itemize-changes", "--stats", "-e", shlex.join(SSH)]
        if transport.protected_args:
            command.append("--protect-args")
        if args.dry_run:
            command.append("--dry-run")
        if args.checksum:
            command.append("--checksum")
        subprocess.run([*command, operand, str(local) + "/"], check=True)
    elif args.command == "resolve":
        remote = remote_path("REMOTE_FOUNDRY_ROOT") / args.run_id
        transport.ssh(f"readlink -f -- {shlex.quote(str(remote))}")
    elif args.command == "run":
        transport.ssh(required("REMOTE_CMD"))
    else:
        # Do not print SSH identity details or remote paths in the diagnostic summary.
        transport.operand("/")
        transport.ssh("command -v rsync >/dev/null && rsync --version >/dev/null")
        print("Noninteractive SSH and remote rsync: OK")
        for name in ("REMOTE_FOLDING_ROOT", "REMOTE_FOUNDRY_ROOT", "REMOTE_CAMPAIGN_ROOT", "REMOTE_LOG_ROOT"):
            if os.environ.get(name):
                path = shlex.quote(str(remote_path(name)))
                transport.ssh(f"test -d {path} && test -r {path}")
                print(f"{name}: readable")
        print("Remote argument mode: " + ("protected" if transport.protected_args else "legacy shell quoting"))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        if isinstance(error, subprocess.CalledProcessError):
            # Avoid echoing private paths/commands into reusable reports.
            print(f"Sync failed: {Path(error.cmd[0]).name} exited with status {error.returncode}.", file=sys.stderr)
            sys.exit(error.returncode if error.returncode > 0 else 1)
        print(f"Sync failed: {error}", file=sys.stderr)
        sys.exit(2)
