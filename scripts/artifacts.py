#!/usr/bin/env python3
"""Resolve recorded project artifacts without changing project state."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys

STAGES = ('01_staging', '02_design', '03_filtering', '04_analysis')
IDENTITY = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


def _relative(value, label):
    if not isinstance(value, str) or not value:
        raise ValueError(f'{label} must be a nonempty relative path')
    path = PurePosixPath(value)
    if path.is_absolute() or '..' in path.parts or not path.parts or str(path) in ('', '.'):
        raise ValueError(f'{label} must be relative and contain no parent traversal')
    return path


def _map(project):
    file = project / 'local/artifact-map.json'
    if not file.exists():
        return [], {}
    data = json.loads(file.read_text())
    if not isinstance(data, dict) or data.get('schema_version') != 1:
        raise ValueError('Unsupported artifact map schema')
    roots = data.get('source_roots')
    moves = data.get('moves')
    if not isinstance(roots, list) or not isinstance(moves, dict):
        raise ValueError('Artifact map requires source_roots and moves')
    aliases = []
    for value in roots:
        if not isinstance(value, str) or not value:
            raise ValueError('Artifact source roots must be absolute paths')
        raw_root = Path(value).expanduser()
        if '..' in raw_root.parts:
            raise ValueError('Artifact source roots must not contain parent traversal')
        root = Path(os.path.normpath(str(raw_root)))
        if not root.is_absolute():
            raise ValueError('Artifact source roots must be absolute paths without parent traversal')
        aliases.append(root)
    if len(set(aliases)) != len(aliases):
        raise ValueError('Artifact source roots must be unique')
    checked = {}
    for source, target in moves.items():
        source_path = _relative(source, 'Artifact move source')
        target_path = _relative(target, 'Artifact move target')
        checked[source_path] = target_path
    return aliases, checked


def _under(path, root):
    try:
        return path.relative_to(root)
    except ValueError:
        return None


def _mapped(relative, moves):
    candidates = [(source, target) for source, target in moves.items()
                  if relative == source or source in relative.parents]
    if not candidates:
        return None
    source, target = max(candidates, key=lambda pair: len(pair[0].parts))
    suffix = relative.relative_to(source)
    return target / suffix


def resolve(project, path):
    """Resolve an artifact to its canonical project path or explicit external path."""
    project = Path(project).expanduser().resolve()
    value = Path(path).expanduser()
    if '..' in value.parts:
        raise ValueError('Artifact path must not contain parent traversal')
    aliases, moves = _map(project)

    if value.is_absolute():
        lexical = Path(os.path.normpath(str(value)))
        canonical = _under(lexical, project)
        if canonical is not None and canonical.parts and canonical.parts[0] == 'campaigns':
            return project / canonical
        matching = [root for root in [project, *aliases]
                    if _under(lexical, root) is not None]
        if not matching:
            return value.resolve()
        alias = max(matching, key=lambda root: len(root.parts))
        relative = _under(lexical, alias)
    else:
        relative = _relative(value.as_posix(), 'Artifact path')
        if relative.parts[0] == 'campaigns':
            return project / Path(*relative.parts)

    mapped = _mapped(relative, moves)
    if mapped is not None:
        destination = project / Path(*mapped.parts)
        if _under(destination, project) is None or _under(destination.resolve(), project) is None:
            raise ValueError('Artifact mapping escapes the project')
        return destination

    candidate = project / Path(*relative.parts)
    if candidate.exists():
        return candidate
    if relative.parts[0] == 'reference':
        legacy = candidate
    else:
        legacy = project / 'reference' / Path(*relative.parts)
    if _under(legacy, project) is None:
        raise ValueError('Artifact path escapes the project')
    return legacy


def locate(project, path, verify=False):
    resolved = resolve(project, path)
    present = resolved.exists()
    digest = None
    if verify and resolved.is_file():
        h = hashlib.sha256()
        with resolved.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(block)
        digest = h.hexdigest()
    return {'path': str(resolved), 'present': present, 'sha256': digest}


def _iter_artifact_entries(project, campaign, stage=None):
    if not IDENTITY.fullmatch(campaign):
        raise ValueError('Use a simple campaign ID')
    project = Path(project).expanduser().resolve()
    campaign_dir = project / 'campaigns' / campaign
    file = campaign_dir / 'artifacts.jsonl'
    if file.is_symlink() or campaign_dir.is_symlink() or (project / 'campaigns').is_symlink():
        raise ValueError('Artifact catalog paths must not be symlinks')
    if not file.is_file():
        raise FileNotFoundError(f'Artifact catalog not found for campaign {campaign}')
    with file.open() as stream:
        for number, raw_line in enumerate(stream, 1):
            line = raw_line.rstrip('\r\n')
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f'Invalid artifact entry on line {number}') from exc
            if not isinstance(record, dict):
                raise ValueError(f'Artifact entry on line {number} must be an object')
            if stage is None or record.get('stage') == stage:
                yield line, record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    locate_parser = commands.add_parser('locate')
    locate_parser.add_argument('path')
    locate_parser.add_argument('--verify', action='store_true')
    locate_parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parent.parent)
    list_parser = commands.add_parser('list')
    list_parser.add_argument('--campaign', required=True)
    list_parser.add_argument('--stage', choices=STAGES)
    list_parser.add_argument('--files', action='store_true')
    list_parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args(argv)
    try:
        if args.command == 'locate':
            print(json.dumps(locate(args.project, args.path, args.verify), separators=(',', ':')))
        elif args.files:
            for line, _ in _iter_artifact_entries(args.project, args.campaign, args.stage):
                print(line)
        else:
            stages = {}
            count = 0
            for _, record in _iter_artifact_entries(args.project, args.campaign, args.stage):
                name = record.get('stage', 'unspecified')
                stages[name] = stages.get(name, 0) + 1
                count += 1
            print(json.dumps({'campaign_id': args.campaign, 'count': count,
                              'stages': stages}, separators=(',', ':'), sort_keys=True))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f'artifacts: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
