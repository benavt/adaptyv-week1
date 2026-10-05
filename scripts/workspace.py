#!/usr/bin/env python3
"""Project front door to the pinned shared workflow implementation."""
import json
import os
from pathlib import Path
import runpy
import sys

project=Path(__file__).resolve().parent.parent
config=project/'local/fleet.json'
private=json.loads(config.read_text()) if config.is_file() else {}
# The private binding selects a runtime without changing shell activation.
selected=private.get('control_python')
if selected:
    interpreter=Path(selected).expanduser()
    if not interpreter.is_absolute() or not interpreter.is_file() or not os.access(interpreter,os.X_OK):
        raise SystemExit('Configured control interpreter is unavailable; inspect Fleet environment binding.')
    if Path(sys.executable).resolve()!=interpreter.resolve() or not sys.flags.isolated:
        os.execv(str(interpreter),[str(interpreter),'-I','-B',str(Path(__file__).resolve()),*sys.argv[1:]])
# Historical manifests are navigation evidence, not reconstructed frozen runs.
if sys.argv[1:2]==['history']:
    import hashlib
    import re
    arguments=sys.argv[2:]
    if len(arguments)>1:
        raise SystemExit('Usage: workspace.py history [CAMPAIGN]')
    selected_campaign=arguments[0] if arguments else None
    if selected_campaign and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',selected_campaign):
        raise SystemExit('Invalid historical campaign ID.')
    base=project/'campaigns'
    sources=[base/selected_campaign/'04_analysis/history/manifest.json'] if selected_campaign else sorted(base.glob('*/04_analysis/history/manifest.json'))
    records=[]
    for source in sources:
        if not source.is_file() or not source.resolve().is_relative_to(base.resolve()):
            raise SystemExit('Historical campaign manifest unavailable.')
        record=json.loads(source.read_text())
        if record.get('schema_version')!=1 or record.get('campaign_id')!=source.parents[2].name:
            raise SystemExit('Historical campaign identity/schema mismatch.')
        item={'campaign_id':record['campaign_id'],'title':record['title'],'summary':record['summary'],
              'manifest_path':source.relative_to(project).as_posix(),
              'manifest_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'home':source.parent.relative_to(project).as_posix(),
              'recorded_run_count':len(record.get('runs',[])),
              'current_output_count':len(record.get('current_outputs',[]))}
        if selected_campaign:
            item.update({key:record.get(key,[]) for key in ['locations','current_outputs','runs','notes','documented_unavailable']})
        records.append(item)
    print(json.dumps({'project_id':'EGFR','historical_root':'campaigns','kind':'historical_navigation','campaigns':records,
          'path_resolver':'python3 scripts/artifacts.py locate ORIGINAL_PATH --verify',
          'policy':'Original scientific IDs and manifest bytes preserved; use the resolver for recorded historical paths. Status is curated local evidence; historical manifests are not immutable run snapshots or live remote monitoring.'},indent=2))
    raise SystemExit(0)
home=os.environ.get('WORKFLOWS_HOME') or private.get('workflow_home')
if not home:raise SystemExit('Set WORKFLOWS_HOME or configure ignored local/fleet.json.')
script=Path(home).expanduser()/'scripts/workflow.py'
if not script.is_file():raise SystemExit('Shared workflow unavailable at configured location.')
os.environ['PROJECT_ROOT']=str(project)
sys.path.insert(0,str(script.parent))
runpy.run_path(str(script),run_name='__main__')
