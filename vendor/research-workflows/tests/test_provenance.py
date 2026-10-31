import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import workflow as w
import provenance as p


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        root=Path(self.temp.name).resolve();self.project=root/'project';self.project.mkdir()
        (self.project/'campaigns').mkdir();(self.project/'submissions').mkdir()
        w.write(self.project/'project.json',{'schema_version':1,'project_id':'TNFA','title':'Test','target':{},'competition':{'challenge_number':2},'default_recipe':{'stages':{s:{} for s in w.STAGES}}})
        w.write(self.project/'workflow.lock.json',{'workflow':{'files_sha256':{'scripts/workflow.py':w.sha(w.HOME/'scripts/workflow.py')}},'guides':{}})
        source=root/'native';source.mkdir();(source/'request').write_text('protocol only')
        settings=root/'settings.json';settings.write_text('{"test_only":true}')
        w.new_campaign(self.project,'test','Test','Protocol test')
        self.run=w.freeze(self.project,'test','test',source,settings)
        self.output=w.run_path(self.project,'test',self.run)/'02_design/native.fasta'
        self.output.write_text('>one\nACDEFGHIKL\n>two\nACDEFGHIKM\n>same\nACDEFGHIKL\n')

    def candidate(self,name='one',parents=()):
        return {'candidate_id':name,'sequence':p.fasta_record(self.output,name),'molecule_class':'protein',
                'origin':{'campaign_id':'test','run_id':self.run,'path':'02_design/native.fasta','sha256':w.sha(self.output),'fasta_id':name},
                'generation':{'method':'synthetic-test-only'},'parents':list(parents)}

    def test_exact_sequence_and_parent_trace(self):
        record=self.candidate();record['sequence']='WRONG'
        with self.assertRaises(ValueError):p.register(self.project,record)
        p.register(self.project,self.candidate());p.register(self.project,self.candidate('two',['one']))
        trace=p.trace(self.project,'two');self.assertEqual(trace['parents'][0]['candidate']['candidate_id'],'one')
        self.assertIn('execution_receipt',trace['candidate']['generation_unknowns'])

    def test_source_hash_changes_refused(self):
        p.register(self.project,self.candidate());self.output.write_text('changed')
        with self.assertRaises(ValueError):p.trace(self.project,'one')

    def test_observation_value_verified_and_selection_frozen(self):
        p.register(self.project,self.candidate());source=self.output.parent.parent/'03_filtering/metrics.json';source.write_text('{"score":0.8}')
        record={'observation_id':'score1','candidate_id':'one','method':'synthetic','version':'1','scope':'test prediction',
                'artifact':{'campaign_id':'test','run_id':self.run,'path':'03_filtering/metrics.json','sha256':w.sha(source)},
                'locator':{'format':'json','pointer':'/score'},'value':0.7}
        with self.assertRaises(ValueError):p.observe(self.project,record)
        record['value']=0.8;p.observe(self.project,record)
        selected=p.select(self.project,['one'],{'method':'test','version':'1','rationale':'test-only'})
        before=(selected/'selection.json').read_bytes()
        record['observation_id']='score2';p.observe(self.project,record)
        self.assertEqual((selected/'selection.json').read_bytes(),before)
        self.assertEqual(w.read(selected/'selection.json')['review_status'],'pending_researcher_review')

    def test_duplicate_sequence_rejected_across_origins(self):
        p.register(self.project,self.candidate());p.register(self.project,self.candidate('same'))
        with self.assertRaises(ValueError):p.select(self.project,['one','same'],{'method':'test','version':'1','rationale':'test'})


    def observation(self,source,locator,value):
        return {'observation_id':'score1','candidate_id':'one','method':'synthetic','version':'1','scope':'test prediction',
                'artifact':{'campaign_id':'test','run_id':self.run,'path':source.relative_to(self.output.parent.parent).as_posix(),'sha256':w.sha(source)},
                'locator':locator,'value':value}

    def test_json_locator_and_native_type_are_strict(self):
        p.register(self.project,self.candidate())
        source=self.output.parent.parent/'03_filtering/metrics.json'
        source.write_text('{"items":[true,1],"a/b":{"~x":0.8}}')
        for pointer in ['items/0','/items/-1','/items/01','/items/2','/items/0/x','/a~2b']:
            with self.subTest(pointer=pointer),self.assertRaises(ValueError):
                p.observe(self.project,self.observation(source,{'format':'json','pointer':pointer},1))
        with self.assertRaises(ValueError):
            p.observe(self.project,self.observation(source,{'format':'json','pointer':'/items/0'},1))
        p.observe(self.project,self.observation(source,{'format':'json','pointer':'/a~1b/~0x'},0.8))
        self.assertEqual(p.trace(self.project,'one')['observations'][0]['value'],0.8)

    def test_csv_locator_rejects_missing_or_ambiguous_data(self):
        p.register(self.project,self.candidate())
        source=self.output.parent.parent/'03_filtering/metrics.csv';source.write_text('name,score\none,0.8\n')
        for row,column in [(1,'score'),(3,'score'),('2','score'),(True,'score'),(2,'missing')]:
            with self.subTest(row=row,column=column),self.assertRaises(ValueError):
                p.observe(self.project,self.observation(source,{'format':'csv','row':row,'column':column},'0.8'))
        source.write_text('score,score\n0.8,0.9\n')
        with self.assertRaises(ValueError):
            p.observe(self.project,self.observation(source,{'format':'csv','row':2,'column':'score'},'0.9'))

    def test_trace_rechecks_native_sequence_and_parent_record(self):
        p.register(self.project,self.candidate());p.register(self.project,self.candidate('two',['one']))
        parent_path=p.candidate_path(self.project,'one');record=w.read(parent_path)
        record['generation']['method']='changed';parent_path.chmod(0o600);w.write(parent_path,record)
        with self.assertRaises(ValueError):p.trace(self.project,'two')
        # A self-consistent edited sequence/hash still must match its native FASTA.
        record['sequence']='ACDEFGHIKM';record['sequence_sha256']=p.hashlib.sha256(record['sequence'].encode()).hexdigest()
        w.write(parent_path,record)
        with self.assertRaises(ValueError):p.trace(self.project,'one')
        with self.assertRaises(ValueError):p.register(self.project,self.candidate('same',['one']))

    def test_parent_sequence_reference_and_candidate_identity_checked(self):
        p.register(self.project,self.candidate());p.register(self.project,self.candidate('two',['one']))
        path=p.candidate_path(self.project,'two');record=w.read(path);path.chmod(0o600)
        record['parents'][0]['sequence_sha256']='0'*64;w.write(path,record)
        with self.assertRaises(ValueError):p.trace(self.project,'two')
        record['parents'][0]['sequence_sha256']=w.read(p.candidate_path(self.project,'one'))['sequence_sha256']
        record['project_id']='OTHER';w.write(path,record)
        with self.assertRaises(ValueError):p.trace(self.project,'two')

    def test_trace_rechecks_observation_value(self):
        p.register(self.project,self.candidate())
        source=self.output.parent.parent/'03_filtering/metrics.json';source.write_text('{"score":0.8}')
        p.observe(self.project,self.observation(source,{'format':'json','pointer':'/score'},0.8))
        path=self.project/'registry/observations/one/score1.json';record=w.read(path)
        record['value']=0.7;path.chmod(0o600);w.write(path,record)
        with self.assertRaises(ValueError):p.trace(self.project,'one')

    def test_ranking_evidence_is_explicit_verified_and_frozen(self):
        p.register(self.project,self.candidate())
        source=self.output.parent.parent/'03_filtering/metrics.json';source.write_text('{"score":0.8}')
        observation=self.observation(source,{'format':'json','pointer':'/score'},0.8);p.observe(self.project,observation)
        ranking={'method':'test','version':'1','rationale':'test-only','evidence':[{'candidate_id':'one','observation_id':'absent'}]}
        with self.assertRaises(ValueError):p.select(self.project,['one'],ranking)
        ranking['evidence'][0]['observation_id']='score1'
        selected=p.select(self.project,['one'],ranking);data=w.read(selected/'selection.json')
        self.assertEqual(data['ranking_evidence'][0]['observation']['value'],0.8)
        self.assertEqual(data['ranking_evidence_unknowns'],[])
        before=(selected/'selection.json').read_bytes();observation['observation_id']='score2';p.observe(self.project,observation)
        self.assertEqual((selected/'selection.json').read_bytes(),before)

    def test_score_candidate_binding_mismatch_and_trace_tamper_refused(self):
        p.register(self.project,self.candidate())
        source=self.output.parent.parent/'03_filtering/metrics.json';source.write_text('{"score":0.8,"sequence":"ACDEFGHIKM","native_id":"one"}')
        record=self.observation(source,{'format':'json','pointer':'/score'},0.8)
        record['candidate_binding']={'kind':'sequence','locator':{'format':'json','pointer':'/sequence'}}
        with self.assertRaises(ValueError):p.observe(self.project,record)
        record['candidate_binding']={'kind':'native_fasta_id','locator':{'format':'json','pointer':'/native_id'}}
        p.observe(self.project,record)
        path=self.project/'registry/observations/one/score1.json';stored=w.read(path)
        stored['candidate_binding']['kind']='sequence';stored['candidate_binding']['locator']['pointer']='/sequence'
        path.chmod(0o600);w.write(path,stored)
        with self.assertRaises(ValueError):p.trace(self.project,'one')

    def test_score_candidate_binding_matching_hash_and_explicit_unknowns(self):
        p.register(self.project,self.candidate())
        candidate=w.read(p.candidate_path(self.project,'one'))
        source=self.output.parent.parent/'03_filtering/metrics.json'
        source.write_text(json.dumps({'score':0.8,'sequence_sha256':candidate['sequence_sha256']}))
        record=self.observation(source,{'format':'json','pointer':'/score'},0.8)
        p.observe(self.project,record)
        self.assertTrue(p.trace(self.project,'one')['observations'][0]['score_candidate_binding_unknowns'])
        record['observation_id']='score2';record['candidate_binding']={'kind':'sequence_sha256','locator':{'format':'json','pointer':'/sequence_sha256'}}
        p.observe(self.project,record)
        observations=p.trace(self.project,'one')['observations']
        self.assertEqual(observations[1]['score_candidate_binding_unknowns'],[])
        self.assertIn('source sequence_sha256 verified',observations[1]['score_candidate_binding_validation'])

    def test_registry_and_submission_symlinks_refused(self):
        other=self.project.parent/'other';other.mkdir()
        (self.project/'registry').symlink_to(other,target_is_directory=True)
        with self.assertRaises(ValueError):p.register(self.project,self.candidate())
        self.assertFalse(any(other.iterdir()))



if __name__=='__main__':unittest.main()
