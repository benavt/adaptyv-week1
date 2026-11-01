import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('workflow',Path(__file__).resolve().parents[1]/'scripts/workflow.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name).resolve();self.project=self.root/'project';self.project.mkdir()
        (self.project/'campaigns').mkdir();self.source=self.root/'native';self.source.mkdir()
        (self.source/'request.json').write_text('{"construct":"test-only"}')
        self.settings=self.root/'settings.json';self.settings.write_text('{"seed":5}')
        w.write(self.project/'project.json',{'schema_version':1,'project_id':'TNFA','title':'Test','target':{},'default_recipe':{'stages':{s:{'method':'explicit'} for s in w.STAGES}}})
        w.write(self.project/'workflow.lock.json',{'workflow':{'files_sha256':{'scripts/workflow.py':w.sha(w.HOME/'scripts/workflow.py')}},'guides':{}})
        w.new_campaign(self.project,'test','Test','Protocol validation, no science')

    def run_new(self,parents=()):
        return w.freeze(self.project,'test','hotspots',self.source,self.settings,parents)

    def test_four_stages_frozen_inputs_and_lineage(self):
        first=self.run_new();p=w.run_path(self.project,'test',first)
        self.assertEqual(sorted(q.name for q in p.iterdir() if q.name in w.STAGES),list(w.STAGES))
        before=(p/'run.json').read_bytes();self.settings.write_text('{"seed":6}')
        second=self.run_new(['test/'+first]);record=w.verify(self.project,'test',second)
        self.assertEqual(record['parents'][0]['run_id'],first)
        self.assertEqual((p/'run.json').read_bytes(),before)
        self.assertRegex(first,r'^\d{4}-\d{2}-\d{2}_hotspots_001$')
        self.assertFalse(record['scientific_inputs_validated'])

    def test_tamper_and_additions_fail(self):
        run=self.run_new();p=w.run_path(self.project,'test',run)/w.STAGES[0]
        p.chmod(0o700);(p/'extra.json').write_text('{}')
        with self.assertRaises(ValueError):w.verify(self.project,'test',run)

    def test_broad_or_managed_source_rejected(self):
        for source in [self.project,self.root,self.project/'campaigns']:
            with self.assertRaises(ValueError):w.freeze(self.project,'test','bad',source,self.settings)

    def test_guide_pin_and_parent_reference_rejected(self):
        lock=w.read(self.project/'workflow.lock.json');lock['workflow']['files_sha256']['scripts/workflow.py']='wrong'
        w.write(self.project/'workflow.lock.json',lock)
        with self.assertRaises(ValueError):self.run_new()
        with self.assertRaises(ValueError):w.new_campaign(self.project,'../bad','Bad','Bad')

    def test_source_symlink_and_credentials_rejected(self):
        (self.source/'linked').symlink_to(self.source/'request.json')
        with self.assertRaises(ValueError):self.run_new()
        (self.source/'linked').unlink();(self.source/'.env.local').write_text('private')
        with self.assertRaises(ValueError):self.run_new()

    def remote(self,root,expected,mode):
        return subprocess.run([sys.executable,'-B','-c',w.REMOTE,str(root),expected,mode],capture_output=True,text=True)

    def test_remote_collision_receipt_and_integrity(self):
        run=self.run_new();p=w.run_path(self.project,'test',run);source=p/w.STAGES[0]
        import shutil
        target=self.root/'remote'/w.STAGES[0];target.parent.mkdir()
        self.assertEqual(self.remote(target,w.sha(source/'snapshot.json'),'preflight').stdout.strip(),'upload_required')
        shutil.copytree(source,target);target.chmod(0o700)
        result=self.remote(target,w.sha(source/'snapshot.json'),'verify');self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual((target.parent/'run.json').read_bytes(),(p/'run.json').read_bytes())
        self.assertNotEqual(self.remote(target,'different','preflight').returncode,0)
        target.chmod(0o700);(target/'extra').write_text('wrong')
        self.assertNotEqual(self.remote(target,w.sha(source/'snapshot.json'),'verify').returncode,0)

    def write_fleet_site(self, remote_workspace='/srv/research/projects/adaptyv-week1', registry=True):
        fleet=self.root/'fleet';(fleet/'local/workspaces').mkdir(parents=True)
        (fleet/'local').mkdir(exist_ok=True)
        w.write(fleet/'local/machines.json',{'machines':{'remote-gpu':{'ssh_alias':'remote-gpu'}}})
        w.write(fleet/'local/workspaces/remote-gpu.json',{'root':'/srv/research'})
        w.write(self.project/'local/fleet.json',{'fleet_home':str(fleet)})
        config=w.read(self.project/'project.json');config['machine_id']='remote-gpu';w.write(self.project/'project.json',config)
        if registry:
            w.write(fleet/'local/projects.json',{'schema_version':1,'projects':{
                'TNFA':{'project_id':'TNFA','machine_id':'remote-gpu','remote_workspace':remote_workspace}}})
        return fleet

    def test_site_uses_registry_slug_and_rejects_invalid_mapping(self):
        fleet=self.write_fleet_site()
        _,_,remote=w.site(self.project)
        self.assertEqual(remote,Path('/srv/research/projects/adaptyv-week1'))
        self.assertEqual(w.project_config(self.project)['project_id'],'TNFA')
        (fleet/'local/projects.json').unlink()
        self.assertEqual(w.site(self.project)[2],Path('/srv/research/projects/TNFA'))
        w.write(fleet/'local/projects.json',{'schema_version':1,'projects':{
            'TNFA':{'project_id':'TNFA','machine_id':'remote-gpu','remote_workspace':'/srv/research/projects/adaptyv-week1'}}})
        registry=fleet/'local/projects.json'
        for entry in [
            {'project_id':'TNFA','machine_id':'remote-gpu','remote_workspace':'/srv/research/elsewhere'},
            {'project_id':'TNFA','machine_id':'other','remote_workspace':'/srv/research/projects/adaptyv-week1'},
        ]:
            with self.subTest(entry=entry):
                w.write(registry,{'schema_version':1,'projects':{'TNFA':entry}})
                with self.assertRaises(ValueError):w.site(self.project)

    def test_output_selection_refuses_inputs_and_overwrite(self):
        run=self.run_new();selected=self.root/'list.txt';selected.write_text('01_staging/payload/request.json\n')
        with patch.object(w,'site',return_value=(self.root,{'ssh_alias':'example'},Path('/example/projects/TNFA'))):
            with self.assertRaises(ValueError):w.transfer(self.project,'test',run,selected=selected)
            selected.write_text('03_filtering/scores.csv\n')
            result=w.transfer(self.project,'test',run,selected=selected);self.assertFalse(result['applied'])
            output=w.run_path(self.project,'test',run)/'03_filtering/scores.csv';output.write_text('keep')
            with self.assertRaises(ValueError):w.transfer(self.project,'test',run,selected=selected)

    def test_rsync_owner_only_permissions(self):
        import shutil
        if not shutil.which('rsync'):self.skipTest('rsync unavailable')
        output=self.root/'copied'
        subprocess.run(['rsync','-rlt','--checksum',str(self.source)+'/',str(output)],check=True,umask=0o077)
        w.protect_outputs(output,['request.json'])
        self.assertEqual(output.stat().st_mode & 0o077,0)
        self.assertEqual((output/'request.json').stat().st_mode & 0o077,0)

    def test_remote_symlink_refused(self):
        link=self.root/'remote-link';link.symlink_to(self.source,target_is_directory=True)
        self.assertNotEqual(self.remote(link,'unknown','preflight').returncode,0)


if __name__=='__main__':unittest.main()
