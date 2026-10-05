#!/usr/bin/env python3
"""Scan the current directory and open the inventory viewer."""

from __future__ import annotations

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
APP_DIR = REPO_ROOT / "desktop_inventory"
PORT = 8876
LOCAL_MODULES = {"server", "scanner"}


def setup_commands(missing_package: str | None = None) -> str:
    lines = [
        "python3 --version",
        f"cd {REPO_ROOT}",
        "python3 inventory.py",
    ]
    if missing_package:
        lines.insert(0, f"python3 -m pip install {missing_package}")
    return (
        "Python 3.10 or newer is required. desktop_inventory uses only the standard library.\n\n"
        + "\n".join(lines)
    )


def fail(message: str, missing_package: str | None = None) -> None:
    print(message, file=sys.stderr)
    print(file=sys.stderr)
    print(setup_commands(missing_package), file=sys.stderr)
    raise SystemExit(1)


def missing_package_name(exc: BaseException) -> str | None:
    if not isinstance(exc, ModuleNotFoundError):
        return None
    name = (exc.name or "").split(".", 1)[0]
    if not name or name in LOCAL_MODULES:
        return None
    return name


def main() -> None:
    if sys.version_info < (3, 10):
        fail(f"Python {sys.version_info.major}.{sys.version_info.minor} is too old.")
    if not APP_DIR.is_dir():
        fail(f"Missing folder: {APP_DIR}")
    sys.path.insert(0, str(APP_DIR))
    try:
        from server import main as serve
    except (ImportError, ModuleNotFoundError) as exc:
        package = missing_package_name(exc)
        if package is None:
            fail(f"Could not import the inventory app from {APP_DIR}: {exc}")
        fail(f"Could not start the inventory server: {exc}", package)

    sys.argv = [sys.argv[0], "--root", os.getcwd(), "--port", str(PORT)]
    serve()


if __name__ == "__main__":
    main()
