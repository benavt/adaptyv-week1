"""Sequence-level lineage and frozen cross-campaign submission selections."""
import csv
import datetime as dt
import hashlib
import json
import re
from pathlib import Path
from zoneinfo import ZoneInfo
import workflow as w


def artifact(project,record):
    path=Path(record['path'])
    if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0] not in w.STAGES:
        raise ValueError('Artifact path must be run-relative within a pipeline stage')
    root=w.run_path(project,record['campaign_id'],record['run_id'])
    w.verify(project,record['campaign_id'],record['run_id'])
    source=root/path
    for item in [source,*source.parents]:
        if item==root.parent:break
        if item.is_symlink():raise ValueError('Provenance artifact is a symlink')
    if not source.is_file() or w.sha(source)!=record['sha256']:
        raise ValueError('Provenance source absent or its hash changed')
    return source


def fasta_record(path,identity):
    records={};name=None;parts=[]
    for line in path.read_text().splitlines()+['>__END__']:
        if line.startswith('>'):
            if name is not None:
                if name in records:raise ValueError('Ambiguous repeated FASTA record identity')
                records[name]=''.join(parts)
            tokens=line[1:].split()
            if not tokens:raise ValueError('Empty FASTA record identity')
            name=tokens[0];parts=[]
        elif name is not None:parts.append(line.strip())
    if identity not in records:raise ValueError('Native FASTA record identity absent')
    return records[identity]


def candidate_path(project,identity):
    w.identity(identity)
    path=project/'registry/candidates'/(identity+'.json')
    managed_path(project,path)
    return path


def managed_path(project,path):
    for item in [path,*path.parents]:
        if item==project.parent:break
        if item.is_symlink():raise ValueError('Provenance registry path is a symlink')


def native_value(source,locator):
    """Resolve a strict RFC 6901 pointer or a numbered CSV data row."""
    if locator.get('format')=='json':
        pointer=locator.get('pointer')
        if not isinstance(pointer,str) or (pointer and not pointer.startswith('/')):
            raise ValueError('JSON pointer must be empty or start with /')
        actual=w.read(source)
        for part in pointer[1:].split('/') if pointer else []:
            if re.search(r'~(?![01])',part):raise ValueError('Malformed JSON pointer escape')
            part=part.replace('~1','/').replace('~0','~')
            if isinstance(actual,list):
                if not re.fullmatch(r'0|[1-9][0-9]*',part):raise ValueError('JSON array index must be a nonnegative canonical integer')
                index=int(part)
                if index>=len(actual):raise ValueError('JSON array index out of range')
                actual=actual[index]
            elif isinstance(actual,dict):
                if part not in actual:raise ValueError('JSON pointer key absent')
                actual=actual[part]
            else:raise ValueError('JSON pointer descends through a scalar')
        return actual
    if locator.get('format')=='csv':
        row=locator.get('row');column=locator.get('column')
        if type(row) is not int or row<2:raise ValueError('CSV source row must be an integer starting at 2 after the header')
        with source.open(newline='') as stream:
            reader=csv.DictReader(stream)
            fields=reader.fieldnames or []
            if len(fields)!=len(set(fields)):raise ValueError('CSV column headers are ambiguous')
            if column not in fields:raise ValueError('CSV source column absent')
            for number,values in enumerate(reader,2):
                if number==row:
                    if None in values or values[column] is None:raise ValueError('CSV source row has inconsistent columns')
                    return values[column]
        raise ValueError('CSV source row out of range')
    raise ValueError('Use an explicit JSON pointer or CSV row/column locator')


def same_value(actual,recorded):
    # Python considers True == 1 and nested values inherit that ambiguity.
    if type(actual) is not type(recorded):return False
    if isinstance(actual,dict):return actual.keys()==recorded.keys() and all(same_value(actual[k],recorded[k]) for k in actual)
    if isinstance(actual,list):return len(actual)==len(recorded) and all(same_value(a,b) for a,b in zip(actual,recorded))
    return actual==recorded



def candidate_binding(source,candidate,record):
    binding=record.get('candidate_binding')
    if binding is None:
        return ('researcher-declared association, not certified',
                ['No source locator links this score to the candidate sequence or native FASTA identity'])
    if not isinstance(binding,dict):raise ValueError('Candidate binding must specify kind and source locator')
    kind=binding.get('kind')
    expected={'sequence':candidate['sequence'],'sequence_sha256':candidate['sequence_sha256'],
              'native_fasta_id':candidate['origin']['fasta_id']}
    if kind not in expected:raise ValueError('Candidate binding kind must be sequence, sequence_sha256 or native_fasta_id')
    if not same_value(native_value(source,binding['locator']),expected[kind]):
        raise ValueError('Source score candidate binding differs from the candidate sequence/identity')
    if kind=='native_fasta_id':
        return ('source native FASTA identity verified; sequence association not independently certified',
                ['Source binding exposes a native FASTA name without the sequence or its hash'])
    return ('source '+kind+' verified; computational result is not experimental certification',[])


def validate_candidate(project,identity):
    candidate=w.read(candidate_path(project,identity))
    if candidate['candidate_id']!=identity or candidate['project_id']!=w.project_config(project)['project_id']:
        raise ValueError('Candidate identity or project differs')
    source=artifact(project,candidate['origin'])
    if fasta_record(source,candidate['origin']['fasta_id'])!=candidate['sequence']:
        raise ValueError('Candidate sequence differs from its native source record')
    if w.sha(w.run_path(project,candidate['origin']['campaign_id'],candidate['origin']['run_id'])/'run.json')!=candidate['origin_run_sha256']:
        raise ValueError('Origin run metadata changed')
    if hashlib.sha256(candidate['sequence'].encode()).hexdigest()!=candidate['sequence_sha256']:
        raise ValueError('Candidate sequence/hash changed')
    return candidate


def register(project,record):
    identity=w.identity(record['candidate_id']);sequence=record['sequence']
    cls=record['molecule_class']
    if cls not in {'protein','nanobody','scfv','fab_kappa','fab_lambda'}:raise ValueError('Unsupported molecule class')
    allowed=set('ACDEFGHIKLMNPQRSTVWY')
    chains=sequence.split(':')
    if any(not chain or set(chain)-allowed for chain in chains):raise ValueError('Supply exact canonical amino-acid sequences')
    if cls.startswith('fab_') and len(chains)!=2:raise ValueError('Fab submission needs VH:VL')
    if not cls.startswith('fab_') and len(chains)!=1:raise ValueError('This molecule class requires one submitted chain')
    source=artifact(project,record['origin'])
    if fasta_record(source,record['origin']['fasta_id'])!=sequence:raise ValueError('Candidate sequence differs from its native source record')
    generation=record.get('generation',{})
    if not generation.get('method'):raise ValueError('Record generation method; unknown fields must remain explicit')
    parents=[]
    parent_ids=record.get('parents',[])
    if identity in parent_ids or len(parent_ids)!=len(set(parent_ids)):raise ValueError('Parents must be distinct existing candidate identities, excluding self')
    for parent in parent_ids:
        original=trace(project,parent)['candidate'];p=candidate_path(project,parent)
        parents.append({'candidate_id':parent,'record_sha256':w.sha(p),'sequence_sha256':original['sequence_sha256']})
    data={**record,'schema_version':1,'project_id':w.project_config(project)['project_id'],
          'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),'parents':parents,
          'origin_run_sha256':w.sha(w.run_path(project,record['origin']['campaign_id'],record['origin']['run_id'])/'run.json'),
          'generation_unknowns':[key for key in ['version','parameters','input_artifacts','execution_receipt'] if not generation.get(key)],
          'recorded_at':dt.datetime.now(dt.timezone.utc).isoformat(),
          'validation':'native FASTA identity/hash verified; generation metadata is recorded, not independently certified'}
    w.write(candidate_path(project,identity),data,exclusive=True,mode=0o400)
    return identity


def observe(project,record):
    candidate=trace(project,record['candidate_id'])['candidate']
    identity=w.identity(record['observation_id']);source=artifact(project,record['artifact'])
    if not record.get('method') or not record.get('version') or not record.get('scope'):
        raise ValueError('Record method, version and scientific scope')
    if not same_value(native_value(source,record['locator']),record['value']):
        raise ValueError('Observation value differs from the source; preserve the native type/value')
    binding_validation,binding_unknowns=candidate_binding(source,candidate,record)
    path=project/'registry/observations'/candidate['candidate_id']/(identity+'.json')
    managed_path(project,path)
    w.write(path,{**record,'schema_version':1,'sequence_sha256':candidate['sequence_sha256'],
                 'candidate_record_sha256':w.sha(candidate_path(project,candidate['candidate_id'])),
                 'score_candidate_binding_validation':binding_validation,
                 'score_candidate_binding_unknowns':binding_unknowns,
                 'recorded_at':dt.datetime.now(dt.timezone.utc).isoformat()},exclusive=True,mode=0o400)
    return identity


def trace(project,identity,seen=None):
    seen=set() if seen is None else seen
    if identity in seen:raise ValueError('Candidate lineage cycle')
    seen=seen|{identity};candidate=validate_candidate(project,identity)
    parents=[]
    for link in candidate['parents']:
        p=candidate_path(project,link['candidate_id'])
        if w.sha(p)!=link['record_sha256']:raise ValueError('Parent candidate record changed')
        parent=trace(project,link['candidate_id'],seen)
        if parent['candidate']['sequence_sha256']!=link['sequence_sha256']:raise ValueError('Parent sequence reference differs')
        parents.append(parent)
    observations=[]
    observation_root=project/'registry/observations'/identity
    managed_path(project,observation_root)
    for p in sorted(observation_root.glob('*.json')):
        managed_path(project,p);record=w.read(p)
        if record['observation_id']!=p.stem:raise ValueError('Observation record identity differs')
        if record['candidate_id']!=identity or record['sequence_sha256']!=candidate['sequence_sha256']:raise ValueError('Observation candidate identity differs')
        if record['candidate_record_sha256']!=w.sha(candidate_path(project,identity)):raise ValueError('Observed candidate record changed')
        source=artifact(project,record['artifact'])
        if not same_value(native_value(source,record['locator']),record['value']):raise ValueError('Observation value differs from its recorded source')
        binding_validation,binding_unknowns=candidate_binding(source,candidate,record)
        if record['score_candidate_binding_validation']!=binding_validation or record['score_candidate_binding_unknowns']!=binding_unknowns:
            raise ValueError('Recorded score candidate binding validation differs')
        observations.append(record)
    return {'candidate':candidate,'parents':parents,'observations':observations}


def select(project,identities,ranking):
    if not identities or len(identities)!=len(set(identities)):raise ValueError('Selection needs distinct candidate identities')
    if not ranking.get('method') or not ranking.get('version') or not ranking.get('rationale'):
        raise ValueError('Record explicit ranking method, version and rationale')
    traces=[trace(project,identity) for identity in identities]
    sequences=[x['candidate']['sequence_sha256'] for x in traces]
    if len(sequences)!=len(set(sequences)):raise ValueError('Duplicate submitted sequences; distinct generation origins are not distinct molecules')
    ranking_evidence=[]
    available={}
    def index_trace(entry):
        candidate_id=entry['candidate']['candidate_id']
        for observation in entry['observations']:
            available[(candidate_id,observation['observation_id'])]=observation
        for parent in entry['parents']:index_trace(parent)
    for entry in traces:index_trace(entry)
    evidence=ranking.get('evidence',[])
    if not isinstance(evidence,list):raise ValueError('Ranking evidence must be a list of candidate/observation references')
    used=set()
    for link in evidence:
        key=(link['candidate_id'],link['observation_id'])
        if key not in available or key in used:raise ValueError('Ranking evidence must reference distinct observations in the selected lineages')
        used.add(key);observation=available[key]
        record_path=project/'registry/observations'/key[0]/(key[1]+'.json')
        ranking_evidence.append({'candidate_id':key[0],'observation_id':key[1],
                                 'record_sha256':w.sha(record_path),'observation':observation})
    date=dt.datetime.now(ZoneInfo('America/New_York')).date().isoformat()
    root=project/'submissions';managed_path(project,root);root.mkdir(mode=0o700,exist_ok=True)
    for n in range(1,10000):
        path=root/f'{date}_selection_{n:03d}'
        try:path.mkdir(mode=0o700);break
        except FileExistsError:continue
    else:raise ValueError('No selection identity available')
    # Snapshot all evidence so later observations do not change this selection.
    data={'schema_version':1,'project_id':w.project_config(project)['project_id'],
          'created_at':dt.datetime.now(dt.timezone.utc).isoformat(),'ranking':ranking,
          'ranking_evidence':ranking_evidence,
          'ranking_evidence_unknowns':[] if ranking_evidence else ['No source observations explicitly linked to the ranking decision'],
          'review_status':'pending_researcher_review','eligibility_status':'not_certified',
          'ordered_candidates':traces,'requirements_snapshot':w.project_config(project)['competition']}
    w.write(path/'selection.json',data,exclusive=True,mode=0o400)
    with (path/'ranked.csv').open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['name','sequence','molecule_class'])
        for entry in traces:
            c=entry['candidate'];writer.writerow([c['candidate_id'],c['sequence'],c['molecule_class']])
    (path/'ranked.csv').chmod(0o400)
    (path/'README.md').write_text('# Selection draft\n\nRanked candidates with frozen origin, parent and scoring evidence. Researcher\nreview, eligibility, track limits and experimental interpretation remain required.\nThis package has not been submitted.\n')
    (path/'README.md').chmod(0o400);path.chmod(0o500)
    return path
