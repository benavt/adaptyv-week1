#!/usr/bin/env python3
"""Campaign/run navigation, frozen inputs and explicit SSH transport. No submission."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from zoneinfo import ZoneInfo

HOME = Path(__file__).resolve().parent.parent
STAGES = ('01_staging', '02_design', '03_filtering', '04_analysis')
REMOTE = r'''
import hashlib,json,os,pathlib,sys
root=pathlib.Path(sys.argv[1]);expected=sys.argv[2];mode=sys.argv[3]
for p in [root,*root.parents]:
 if p.is_symlink():raise RuntimeError('Symlink storage path refused')
marker=root/'snapshot.json'
if not marker.exists():
 if mode=='preflight' and (not root.exists() or not any(root.iterdir())):
  print('upload_required');sys.exit(0)
 raise RuntimeError('Snapshot absent from a nonempty input path')
if hashlib.sha256(marker.read_bytes()).hexdigest()!=expected:raise RuntimeError('Different immutable snapshot: choose a new run')
m=json.loads(marker.read_text());actual={}
for p in root.rglob('*'):
 if p.is_symlink():raise RuntimeError('Input symlink refused')
 if p.is_file() and p!=marker:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1048576),b''):h.update(b)
  actual[p.relative_to(root).as_posix()]=h.hexdigest()
if mode=='preflight':
 print('already_verified' if actual==m['files_sha256'] else 'upload_required');sys.exit(0)
if actual!=m['files_sha256']:raise RuntimeError('Input hashes differ or contain unexpected files')
record=root.parent/'run.json';body=marker.read_bytes()
if record.is_symlink():raise RuntimeError('Receipt symlink refused')
if record.exists():
 if record.read_bytes()!=body:raise RuntimeError('Existing run receipt differs')
else:
 fd=os.open(record,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o400)
 with os.fdopen(fd,'wb') as f:f.write(body)
for name in ['02_design','03_filtering','04_analysis','receipts']:
 p=root.parent/name
 if p.is_symlink():raise RuntimeError('Stage symlink refused')
 p.mkdir(mode=0o700,exist_ok=True)
for p in root.rglob('*'):p.chmod(0o500 if p.is_dir() else 0o400)
root.chmod(0o500)
print('verified')
'''


def read(path):
    return json.loads(path.read_text())


def write(path, data, exclusive=False, mode=0o600):
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    with path.open('x' if exclusive else 'w') as stream:
        json.dump(data, stream, indent=2); stream.write('\n')
    path.chmod(mode)


def identity(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError('Use a simple project/campaign/run ID')
    return value


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): h.update(block)
    return h.hexdigest()


def tree(path):
    result = {}
    for item in sorted(path.rglob('*')):
        if item.is_symlink(): raise ValueError('Materialize scientific input symlinks before freezing')
        if item.is_file(): result[item.relative_to(path).as_posix()] = sha(item)
    return result


def project_config(project):
    data = read(project/'project.json')
    if data.get('schema_version') != 1: raise ValueError('Unsupported project schema')
    identity(data['project_id'])
    return data


def campaign_config(project, campaign):
    identity(campaign)
    data = read(project/'campaigns'/campaign/'campaign.json')
    if data['campaign_id'] != campaign: raise ValueError('Campaign identity mismatch')
    return data


def run_path(project, campaign, run):
    identity(campaign); identity(run)
    path=project/'campaigns'/campaign/'runs'/run
    for item in [path,*path.parents]:
        if item==project.parent:break
        if item.is_symlink():raise ValueError('Managed run path is a symlink')
    return path


def lock_check(project):
    lock = read(project/'workflow.lock.json')
    if lock['workflow']['files_sha256'] != {n:sha(HOME/n) for n in lock['workflow']['files_sha256']}:
        raise ValueError('Shared workflow changed: deliberately update the project pin before new runs')
    for guide in lock.get('guides', {}).values():
        path = Path(guide['path'])
        if path.is_absolute() or '..' in path.parts: raise ValueError('Invalid project-relative guide')
        if tree(project/path) != guide['files_sha256']:
            raise ValueError('Guide reference changed: deliberately update its project pin')
    return lock


def new_campaign(project, campaign, title, objective):
    identity(campaign); data = project_config(project)
    recipe = data['default_recipe']
    path = project/'campaigns'/campaign
    if (project/'campaigns').is_symlink():raise ValueError('Campaign root is a symlink')
    (project/'campaigns').mkdir(mode=0o700,exist_ok=True)
    path.mkdir(mode=0o700)
    write(path/'campaign.json', {'schema_version':1,'campaign_id':campaign,
          'title':title,'objective':objective,'recipe':recipe,
          'created_at':dt.datetime.now(dt.timezone.utc).isoformat()}, exclusive=True)
    (path/'runs').mkdir(mode=0o700)
    overview(project)
    return path


def freeze(project, campaign, label, source, settings, parents=()):
    data = project_config(project); campaign_data = campaign_config(project, campaign)
    if not re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', label): raise ValueError('Use a lowercase hyphen-separated label')
    if source.is_symlink(): raise ValueError('Source directory is a symlink')
    source = source.resolve(); project = project.resolve()
    if not source.is_dir() or source == project or source in project.parents:
        raise ValueError('Select an explicit prepared input directory')
    if project in source.parents and source.relative_to(project).parts[0] != 'drafts':
        raise ValueError('Prepare local inputs under drafts/, outside managed campaigns and reference trees')
    for item in source.rglob('*'):
        relative = item.relative_to(source)
        if any(p in {'.git','.ssh','local'} or p.startswith('.env') for p in relative.parts) or item.suffix in {'.key','.pem'} or item.name.startswith('id_'):
            raise ValueError('Private/config artifact found among inputs')
    source_hashes = tree(source)
    if not source_hashes: raise ValueError('Prepared input directory is empty')
    scientific_settings = read(settings)
    if not isinstance(scientific_settings, dict) or not scientific_settings: raise ValueError('Supply explicit scientific settings as a nonempty JSON object')
    lock = lock_check(project); links = []
    for parent in parents:
        parts = parent.split('/')
        if len(parts) != 2: raise ValueError('Parent reference is CAMPAIGN/RUN')
        p = run_path(project, *parts)/'run.json'
        parent_record = read(p)
        if parent_record['project_id'] != data['project_id']: raise ValueError('Parent belongs to another project')
        links.append({'campaign_id':parts[0],'run_id':parts[1],'manifest_sha256':sha(p)})
    date = dt.datetime.now(ZoneInfo('America/New_York')).date().isoformat()
    root = project/'campaigns'/campaign/'runs';root.mkdir(exist_ok=True)
    run_path(project,campaign,'validation')
    for sequence in range(1,10000):
        run = f'{date}_{label}_{sequence:03d}'; path = root/run
        try: path.mkdir(mode=0o700); break
        except FileExistsError: continue
    else: raise ValueError('No available run ID')
    for stage in STAGES: (path/stage).mkdir(mode=0o700)
    (path/'receipts').mkdir(mode=0o700)
    inputs = path/STAGES[0]
    shutil.copytree(source, inputs/'payload')
    if tree(inputs/'payload') != source_hashes or tree(source) != source_hashes:
        raise ValueError('Inputs changed during snapshot; incomplete directory retained for review')
    write(inputs/'settings.json',scientific_settings)
    write(inputs/'campaign.json',campaign_data)
    write(inputs/'project.json',data)
    write(inputs/'workflow.lock.json',lock)
    manifest = {'schema_version':1,'project_id':data['project_id'],'campaign_id':campaign,
                'run_id':run,'created_at':dt.datetime.now(dt.timezone.utc).isoformat(),
                'parents':links,'recipe':campaign_data['recipe'],'settings':scientific_settings,
                'files_sha256':tree(inputs),'scientific_inputs_validated':False}
    write(inputs/'snapshot.json',manifest,exclusive=True,mode=0o400)
    write(path/'run.json',manifest,exclusive=True,mode=0o400)
    write(path/'state.json',{'observed_at':manifest['created_at'],'stages':{
          stage:('frozen_requires_native_validation' if stage==STAGES[0] else 'not_started') for stage in STAGES}})
    for item in inputs.rglob('*'): item.chmod(0o500 if item.is_dir() else 0o400)
    inputs.chmod(0o500)
    overview(project)
    return run


def verify(project,campaign,run):
    path=run_path(project,campaign,run);record=read(path/'run.json');inputs=path/STAGES[0]
    if (record['project_id'],record['campaign_id'],record['run_id']) != (project_config(project)['project_id'],campaign,run):
        raise ValueError('Run identity differs from its location')
    actual=tree(inputs);actual.pop('snapshot.json',None)
    if actual!=record['files_sha256'] or (inputs/'snapshot.json').read_bytes()!=(path/'run.json').read_bytes():
        raise ValueError('Frozen inputs changed; restore them or create a new run')
    return record


def overview(project):
    data=project_config(project)
    lines=['# '+data['title'],'','Generated from campaign/run records; status is recorded evidence, not live monitoring.','',
           '| Campaign | Run | Staging | Design | Filtering | Analysis |',
           '| --- | --- | --- | --- | --- | --- |']
    for path in sorted((project/'campaigns').glob('*/campaign.json')):
        campaign=read(path); name=campaign['campaign_id'];runs=sorted((path.parent/'runs').glob('*/run.json'))
        if not runs:lines.append(f'| [{name}](campaigns/{name}/campaign.json) | No runs | — | — | — | — |')
        for record in runs:
            run=read(record)['run_id']; state=record.parent/'state.json'
            stages=read(state).get('stages',{}) if state.exists() else {}
            cells=[str(stages.get(s,'unknown')).replace('|','/') for s in STAGES]
            lines.append(f'| {name} | [{run}](campaigns/{name}/runs/{run}/run.json) | '+' | '.join(cells)+' |')
    if len(lines)==6:lines+=['','No campaign has been named yet. Create one when its objective is agreed.']
    (project/'overview.md').write_text('\n'.join(lines)+'\n')


def site(project):
    private=read(project/'local/fleet.json')
    fleet=Path(private['fleet_home'])
    data=project_config(project)
    machines=read(fleet/'local/machines.json')['machines'];machine=machines[data['machine_id']]
    catalog=read(fleet/'local/workspaces'/(data['machine_id']+'.json'))
    catalog_projects=Path(catalog['root'])/'projects'
    registry=fleet/'local/projects.json'
    if registry.exists():
        manifest=read(registry)
        if manifest.get('schema_version') != 1:
            raise ValueError('Unsupported Fleet projects registry schema')
        projects=manifest.get('projects')
        entry=projects.get(data['project_id']) if isinstance(projects,dict) else None
        if not isinstance(entry,dict) or entry.get('project_id') != data['project_id']:
            raise ValueError('Fleet projects registry has no matching project')
        if entry.get('machine_id') != data['machine_id']:
            raise ValueError('Fleet project machine differs from the project configuration')
        workspace=entry.get('remote_workspace')
        if not isinstance(workspace,str) or not workspace:
            raise ValueError('Fleet project has no remote workspace')
        remote_root=Path(workspace)
        if not remote_root.is_absolute() or '..' in remote_root.parts:
            raise ValueError('Fleet remote workspace must be an absolute path without parent traversal')
        catalog_projects=Path(os.path.normpath(str(catalog_projects)))
        remote_root=Path(os.path.normpath(str(remote_root)))
        try:
            relative=remote_root.relative_to(catalog_projects)
        except ValueError as exc:
            raise ValueError('Fleet remote workspace must be beneath the catalog projects root') from exc
        if not relative.parts:
            raise ValueError('Fleet remote workspace must name a project folder beneath the catalog root')
    else:
        remote_root=catalog_projects/data['project_id']
    return fleet,machine,remote_root


def ssh(machine):
    return ['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15',
            '-o','ServerAliveInterval=10','-o','ServerAliveCountMax=2',machine['ssh_alias']]


def remote_check(machine,inputs,expected,mode):
    p=subprocess.run(ssh(machine)+[shlex.join(['python3','-c',REMOTE,str(inputs),expected,mode])],
                     text=True,capture_output=True,timeout=60)
    if p.returncode:raise ValueError('Remote input check refused: '+p.stderr[-300:])
    return p.stdout.strip()


def protect_outputs(root,names):
    for name in names:
        path=root/name
        for item in [path,*path.parents]:
            if item==root.parent:break
            if item.is_symlink():raise ValueError('Collected path is a symlink')
        path.chmod(0o600)
        for parent in path.parents:
            if parent==root.parent:break
            parent.chmod(0o700)


def transfer(project,campaign,run,apply=False,selected=None):
    verify(project,campaign,run);_,machine,root=site(project)
    local=run_path(project,campaign,run);remote=root/'campaigns'/campaign/'runs'/run
    prefix=['rsync','-rlt','--checksum','-e',shlex.join(ssh(machine)[:-1])]
    alias=machine['ssh_alias']
    if selected is None:
        inputs=local/STAGES[0]; target=remote/STAGES[0]
        command=prefix+['--rsync-path','umask 077; mkdir -p '+shlex.quote(str(target))+' && rsync',
                        str(inputs)+'/',alias+':'+str(target)+'/']
        if apply:
            match=remote_check(machine,target,sha(inputs/'snapshot.json'),'preflight')
            if match!='already_verified':subprocess.run(command,check=True,timeout=600,umask=0o077)
            remote_check(machine,target,sha(inputs/'snapshot.json'),'verify')
            stamp=dt.datetime.now(dt.timezone.utc);date=stamp.astimezone(ZoneInfo('America/New_York')).date().isoformat()
            write(local/'receipts'/(date+'_upload_'+stamp.strftime('%H-%M-%S-%fZ')+'.json'),
                  {'kind':'verified_input_upload','observed_at':dt.datetime.now(dt.timezone.utc).isoformat(),
                   'remote_root':str(remote),'snapshot_sha256':sha(inputs/'snapshot.json')},exclusive=True)
    else:
        names=[s.strip() for s in selected.read_text().splitlines() if s.strip()]
        for name in names:
            p=Path(name)
            if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0] not in STAGES[1:]:
                raise ValueError('Select explicit relative output files within design/filtering/analysis')
            if len(p.parts)<2 or name.endswith('/'):raise ValueError('Select files, not entire stage directories')
        if not names:raise ValueError('Selected-file list is empty')
        for name in names:
            for item in [local/name,*(local/name).parents]:
                if item==local.parent:break
                if item.is_symlink():raise ValueError('Local collection path is a symlink')
            if (local/name).exists():raise ValueError('Local output already exists; verify/adopt it deliberately before collecting')
        command=prefix+['--files-from',str(selected.resolve()),alias+':'+str(remote)+'/',str(local)+'/']
        if apply:
            # Remote lstat/hash protects against symlinks and checks every selected file.
            code="import pathlib,sys,json,hashlib; r=pathlib.Path(sys.argv[1]); out={}; "
            code+="\nfor n in json.loads(sys.argv[2]):\n p=r/n\n if any(x.is_symlink() for x in [p,*p.parents]) or not p.is_file():raise RuntimeError('Selected output is absent or a symlink')\n h=hashlib.sha256()\n with p.open('rb') as f:\n  for b in iter(lambda:f.read(1048576),b''):h.update(b)\n out[n]=h.hexdigest()\nprint(json.dumps(out))"
            result=subprocess.run(ssh(machine)+[shlex.join(['python3','-c',code,str(remote),json.dumps(names)])],text=True,capture_output=True,timeout=60)
            if result.returncode:raise ValueError('Selected remote files could not be verified: '+result.stderr[-250:])
            expected=json.loads(result.stdout)
            subprocess.run(command,check=True,timeout=600,umask=0o077)
            protect_outputs(local,names)
            if any(not (local/n).is_file() or sha(local/n)!=expected[n] for n in names):raise ValueError('Collected output hashes differ; no verified receipt retained')
            stamp=dt.datetime.now(dt.timezone.utc);date=stamp.astimezone(ZoneInfo('America/New_York')).date().isoformat()
            write(local/'receipts'/(date+'_collection_'+stamp.strftime('%H-%M-%S-%fZ')+'.json'),
                  {'kind':'verified_selected_collection','observed_at':dt.datetime.now(dt.timezone.utc).isoformat(),
                   'remote_root':str(remote),'files_sha256':expected},exclusive=True)
    return {'applied':apply,'command':shlex.join(command),'remote_run':str(remote)}


def show(project,campaign=None,run=None,stage=None):
    data=project_config(project)
    out={'project_id':data['project_id'],'title':data['title'],'target':data['target'],
         'recipe':data['default_recipe'],'campaigns':[p.parent.name for p in sorted((project/'campaigns').glob('*/campaign.json'))]}
    if campaign:out['campaign']=campaign_config(project,campaign)
    if run:
        if not campaign:raise ValueError('Select a campaign before a run')
        out['run']=verify(project,campaign,run)
    if stage:
        out={'project_id':data['project_id'],'stage':stage,'plan':out.get('campaign',{}).get('recipe',data['default_recipe'])['stages'][stage],
             'instructions':str(HOME/'stages'/(stage+'.md'))}
        for name,route in data.get('guide_references',{}).items():
            if stage in route['stages']:out[name+'_guide']=str(project/route['path']/'README.md')
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,default=Path(os.environ.get('PROJECT_ROOT',os.getcwd())))
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('show');p.add_argument('campaign',nargs='?');p.add_argument('run',nargs='?');p.add_argument('--stage',choices=STAGES)
    sub.add_parser('overview')
    p=sub.add_parser('campaign');p.add_argument('id');p.add_argument('--title',required=True);p.add_argument('--objective',required=True)
    p=sub.add_parser('new');p.add_argument('campaign');p.add_argument('label');p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--settings',type=Path,required=True);p.add_argument('--parent-run',action='append',default=[])
    for name in ['verify','upload','collect']:
        p=sub.add_parser(name);p.add_argument('campaign');p.add_argument('run')
        if name!='verify':p.add_argument('--apply',action='store_true')
        if name=='collect':p.add_argument('--files-from',type=Path,required=True)
    p=sub.add_parser('candidate');p.add_argument('action',choices=['register','observe','trace']);p.add_argument('--record',type=Path);p.add_argument('--id')
    p=sub.add_parser('selection');p.add_argument('--ids-file',type=Path,required=True);p.add_argument('--ranking-record',type=Path,required=True)
    p=sub.add_parser('standard');p.add_argument('name',choices=['figures','presentations','molecular-views'])
    args=parser.parse_args();project=args.project.resolve();project_config(project)
    if args.command in {'candidate','selection'}:
        import provenance
        if args.command=='selection':
            print(provenance.select(project,[s.strip() for s in args.ids_file.read_text().splitlines() if s.strip()],read(args.ranking_record)))
        elif args.action=='trace':
            if not args.id:parser.error('candidate trace requires --id')
            print(json.dumps(provenance.trace(project,args.id),indent=2))
        else:
            if not args.record:parser.error('candidate register/observe requires --record')
            print(getattr(provenance,'register' if args.action=='register' else 'observe')(project,read(args.record)))
    elif args.command=='show':print(json.dumps(show(project,args.campaign,args.run,args.stage),indent=2))
    elif args.command=='overview':overview(project);print(project/'overview.md')
    elif args.command=='campaign':print(new_campaign(project,args.id,args.title,args.objective))
    elif args.command=='new':print(freeze(project,args.campaign,args.label,args.input_dir,args.settings,args.parent_run))
    elif args.command=='verify':verify(project,args.campaign,args.run);print('Frozen inputs verified; native scientific validation is separate.')
    elif args.command=='standard':print((HOME/'standards'/(args.name+'.md')).read_text())
    else:print(json.dumps(transfer(project,args.campaign,args.run,args.apply,getattr(args,'files_from',None)),indent=2))


if __name__=='__main__':
    try: main()
    except (ValueError,KeyError,OSError,subprocess.SubprocessError) as error:
        print('Workflow: '+str(error),file=sys.stderr);sys.exit(1)
