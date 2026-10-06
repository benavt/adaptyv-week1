#!/usr/bin/env python3
"""Publish a small, source-checked presentation copy of the saved test panel.

Run with --source-root (scientific checkout) and --repo-root (task worktree).
Original packages, frozen selections, and prediction settings are never edited.
"""
import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

PANEL = Path('candidates/2026-10-06_readme-panel_001')
OBS = Path('candidates/2026-10-06_analysis/2026-10-06_stats_006/observations.csv')
STATS_MANIFEST = OBS.with_name('manifest.json')
SELECTION = Path('candidates/2026-10-05_Test_Submission_20_sequences.provenance.json')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')

def write_csv(path, rows, fields=None):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def number(value, digits):
    return f'{float(value):.{digits}f}' if value else '—'

def source_model(root, row):
    p = root / 'candidates' / row['package']
    d, species = row['design_id'], row['species']
    if 'ph-redesigns' in row['package']:
        rows = list(csv.DictReader((p / 'all_statistics.csv').open()))
        item = next(x for x in rows if x['design_id'] == d and x['species'] == species)
        model = Path(item['native_model'])
    elif 'ligandai' in row['package']:
        model, = (p / 'models' / d / species).glob('*.cif')
    elif 'alphaprotein' in row['package']:
        model, = (p / 'models').glob('*' + d + '/' + species + '/*.cif')
    elif 'three-candidates' in row['package']:
        model = p / 'models' / f'{d}_{species}_model_0.cif'
    elif 'd46-130' in row['package']:
        model = p / 'models' / f'fusion_{species}.cif'
    else:
        raise ValueError('Unsupported source package: ' + row['package'])
    if row['native_model_sha256']:
        assert sha(model) == row['native_model_sha256'], f'Model hash mismatch: {model}'
    return model

def display_name(name):
    if name.startswith('Ligand_AI_b6a74871_Design_'):
        return 'LigandAI ' + name.rsplit('_', 1)[1]
    if name.startswith('DD1_Human_EGFR_APNovo_'):
        return 'APNovo ' + name.rsplit('_', 1)[1]
    if name == 'rank_05_D46_130_model_1_chain_A':
        return 'D46-130 fusion'
    return name.replace('_phswitch_', ' / ')

def preview_cell(row, prefix):
    s = row['species']
    return (f'<a href="{prefix}{s}.png"><img src="{prefix}{s}.png" width="180" alt="{s} EGFR complex" /></a>'
            f'<br>ipSAE **{number(row["ipSAE"], 3)}** · Kd **{number(row["kd"], 2)} nM**'
            f'<br>pKa {number(row["H37"], 2)} / {number(row["H100"], 2)}'
            f'<br>[PNG]({prefix}{s}.png) · [PSE]({prefix}{s}.pse)')

def main():
    args = argparse.ArgumentParser()
    args.add_argument('--source-root', type=Path, required=True)
    args.add_argument('--repo-root', type=Path, required=True)
    opts = args.parse_args()
    src, repo = opts.source_root.resolve(), opts.repo_root.resolve()
    target = repo / PANEL
    if target.exists():
        raise FileExistsError('Presentation snapshot already exists; use a successor package.')
    selection = json.loads((src / SELECTION).read_text())
    selected = selection['selected']
    assert len(selected) == 20 and selection['submitted'] is False
    assert len({x['sequence_sha256'] for x in selected}) == 20
    stats = json.loads((src / STATS_MANIFEST).read_text())
    assert sha(src / OBS) == stats['outputs']['observations.csv']
    observations = list(csv.DictReader((src / OBS).open()))
    # Validate original source FASTAs and score tables against the completed snapshot.
    for item in selection['sources']:
        p = src / item['path']
        if p.exists():
            assert sha(p) == item['sha256'], f'Selection source mismatch: {p}'
    target.mkdir(parents=True)
    records, panel_rows, submission_rows, files = [], [], [], []
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    parser = MMCIFParser(QUIET=True)
    for idx, chosen in enumerate(selected, 1):
        assert hashlib.sha256(chosen['sequence'].encode()).hexdigest() == chosen['sequence_sha256']
        rows = [r for r in observations if r['sequence'] == chosen['sequence']
                and r['package'] == Path(chosen['source']).parts[1]]
        assert len(rows) == 2 and {r['species'] for r in rows} == {'Human', 'Mouse'}
        rows.sort(key=lambda r: r['species'])
        folder = target / chosen['name']
        folder.mkdir()
        (folder / 'sequence.fasta').write_text('>' + chosen['name'] + '\n' + chosen['sequence'] + '\n')
        model_sources = []
        for row in rows:
            table = src / 'candidates' / row['source_table']
            key = row['source_table']
            assert sha(table) == stats['sources'][key], f'Score table hash mismatch: {key}'
            model = source_model(src, row)
            structure = parser.get_structure('complex', model)[0]
            protein = {c.id: ''.join(seq1(r.resname) for r in c if r.id[0] == ' ')
                       for c in structure}
            assert protein['A'] == chosen['sequence'], f'Binder mismatch: {model}'
            receptor = protein['B']
            assert len(receptor) == 192 and receptor[36] == 'H' and receptor[99] == 'H'
            shutil.copyfile(model, folder / (row['species'] + '.cif'))
            model_sources.append({'species': row['species'], 'source_model': str(model.relative_to(src)),
                                  'source_sha256': sha(model), 'binder_chain': 'A', 'egfr_chain': 'B',
                                  'egfr_sequence': receptor,
                                  'egfr_sequence_sha256': hashlib.sha256(receptor.encode()).hexdigest()})
            panel_rows.append(row)
        write_csv(folder / 'metrics.csv', rows)
        detail = {'panel_position': idx, 'candidate_id': chosen['name'], 'score_design_id': rows[0]['design_id'],
                  'sequence_sha256': chosen['sequence_sha256'], 'source_package': rows[0]['package'],
                  'source_score_tables': [{'path': r['source_table'], 'row': r['source_row'],
                                          'sha256': stats['sources'][r['source_table']]} for r in rows],
                  'models': model_sources, 'status': 'placeholder_test_panel_not_submitted',
                  'affinity_method': 'DeltaForge sequence fold-and-score' if 'd46-130' in rows[0]['package']
                                     else 'DeltaForge scoring of supplied structure',
                  'pka_scopes': {r['species']: r['pka_scope'] for r in rows},
                  'source_selection_manifest': str(SELECTION), 'source_selection_sha256': sha(src / SELECTION),
                  'source_statistics_snapshot': str(STATS_MANIFEST), 'source_statistics_manifest_sha256': sha(src / STATS_MANIFEST)}
        dump(folder / 'provenance.json', detail)
        records.append(detail)
        subtitle = f'# {chosen["name"]}\n\nPlaceholder test-panel position **{idx:02d} of 20**; not an affinity rank or a confirmed submission.\n\n'
        subtitle += '[Back to the panel](../README.md) · [Exact sequence](sequence.fasta) · [Metrics CSV](metrics.csv) · [Source evidence](provenance.json)\n\n'
        subtitle += '| Species | Predicted complex | ipSAE | EGFR H37 pKa | EGFR H100 pKa | DeltaForge Kd (nM) | Files |\n| --- | --- | ---: | ---: | ---: | ---: | --- |\n'
        for r in rows:
            s = r['species']
            subtitle += (f'| {s} | <img src="{s}.png" width="280" alt="{s} EGFR complex" /> | {number(r["ipSAE"],3)} | {number(r["H37"],3)} | '
                         f'{number(r["H100"],3)} | {number(r["kd"],3)} | [PSE]({s}.pse) · [PNG]({s}.png) · [CIF]({s}.cif) |\n')
        subtitle += '\nCyan: binder chain A. Grey: EGFR D3 chain B. Green sticks: EGFR H37/H100 (historical full positions 346/409). Magenta marks the DD2 extension in the D46-130 fusion. Native residue numbering and all native atoms are retained. PSE coordinates are rigidly aligned on EGFR for display; the paired native CIFs are unchanged and were used for scoring.\n\n'
        subtitle += 'These are computational predictions. The PNG/PSE view is a presentation of the scored native model, not a pH-conditioned structure.\n\n'
        subtitle += f'Affinity workflow: **{detail["affinity_method"]}**. Kd is not an experimental measurement or a prediction at a specified pH.\n\n'
        for r in rows:
            subtitle += f'- {r["species"]} pKa context: {r["pka_scope"]}. H37 status: `{r["H37_status"]}`; H100 status: `{r["H100_status"]}`.\n'
            for site in ('H37','H100'):
                for condition in ('acid_deprotonated','acid_protonated'):
                    if r[site + '_' + condition]:
                        subtitle += f'  Conditional {site}, {condition.replace("_"," ")}: {r[site+"_"+condition]}. This assumption is kept separate from the displayed unforced value.\n'
        (folder / 'README.md').write_text(subtitle)
        submission_rows.append({'name': chosen['submission_name'], 'sequence': chosen['sequence'], 'molecule_class': chosen['molecule_class']})
    write_csv(target / 'metrics.csv', panel_rows)
    write_csv(target / 'test-panel.csv', submission_rows)
    dump(target / 'selection.json', {'status': 'placeholder_test_panel_not_submitted', 'selection_method': selection['method'],
                                     'source_manifest': str(SELECTION), 'source_manifest_sha256': sha(src / SELECTION),
                                     'original_csv_sha256': selection['output_sha256'],
                                     'csv_note': 'test-panel.csv reconstructed from the saved selected records; original CSV unavailable locally; no claim of original CSV byte identity.',
                                     'selected': selected})
    dump(target / 'manifest.json', {'schema_version':1, 'created_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
                                   'purpose':'GitHub presentation copy; original scientific records unchanged',
                                   'predecessors':[{'path':str(SELECTION),'sha256':sha(src/SELECTION)},
                                                   {'path':str(STATS_MANIFEST),'sha256':sha(src/STATS_MANIFEST)}],
                                   'panel_count':20, 'species_observations':40, 'submitted':False,
                                   'records':records, 'outputs':{str(p.relative_to(target)):sha(p) for p in target.rglob('*') if p.is_file()}})
    intro = ('## Placeholder panel of 20 designs\n\n'
             'The saved **Test Submission** panel is a random sample of 20 exact, deduplicated sequences. '
             '**It has not been submitted and is not the final scientific selection.** Order matches the saved test panel. '
             'Replace it with a reviewed successor panel when the final set is frozen.\n\n'
             'Each species cell shows **ipSAE · predicted DeltaForge Kd in nM**, then **EGFR H37 / H100 pKa**. '
             'H37/H100 are EGFR D3 sites (historical full positions 346/409), not binder histidines. '
             'pKa scopes differ by campaign; candidate pages preserve context, status, and conditional results. '
             'These models and scores do not establish pH-selective binding.\n\n')
    table = '| # | Design and evidence | Human EGFR D3 | Mouse EGFR D3 |\n| ---: | --- | --- | --- |\n'
    for idx, (chosen, record) in enumerate(zip(selected, records), 1):
        rows = [r for r in panel_rows if r['package'] == record['source_package'] and r['design_id'] == record['score_design_id']]
        prefix = PANEL.as_posix() + '/' + chosen['name'] + '/'
        cells = [preview_cell(r, prefix) for r in rows]
        table += f'| {idx:02d} | [{display_name(chosen["name"])}]({prefix}) | {cells[0]} | {cells[1]} |\n'
    local_table = table.replace(PANEL.as_posix() + '/', '')
    (target / 'README.md').write_text('# EGFR test panel — 20 designs\n\n' + intro +
                                     '[Panel CSV](test-panel.csv) · [All metrics](metrics.csv) · [Selection](selection.json) · [Manifest](manifest.json)\n\n' + local_table)
    prior_readme = (repo / 'README.md').read_text()
    prior_readme = prior_readme.replace('(overview.md)', '(../overview.md)').replace('(docs/WORKSPACE.md)', '(WORKSPACE.md)').replace('(docs/CHALLENGE.md)', '(CHALLENGE.md)').replace('(docs/HISTORY.md)', '(HISTORY.md)')
    (repo / 'docs/REPOSITORY.md').write_text(prior_readme)
    root_readme = ('# pH-sensitive EGFR binder design\n\n'
                   'Computational design of binders to Domain III of **human and mouse EGFR**, with a proposed proton-linked DD1/DD2 masking mechanism.\n\n'
                   '![Proposed pH-dependent DD1/DD2 masking mechanism](publication/figures/2026-10-06_egfr-ph-switch-schematic.svg)\n\n'
                   '*Design hypothesis:* at pH 7.4 a tethered DD2 mask competes with EGFR for DD1; at pH 6.5 protonation could weaken the mask and favor DD1–EGFR contacts. '
                   'The cartoon illustrates intended states, not experimentally established structures or affinities.\n\n'
                   '**[Read the manuscript draft](publication/2026-10-06_manuscript/README.md)** · '
                   '[Methods](publication/2026-10-06_manuscript/03-methods.md) · '
                   '[References and original PDFs](publication/2026-10-06_reference-library.md) · '
                   '[Repository workflow](docs/REPOSITORY.md)\n\n' + intro +
                   f'[Panel directory]({PANEL.as_posix()}/) · [Sequence CSV]({PANEL.as_posix()}/test-panel.csv) · '
                   f'[Full metrics]({PANEL.as_posix()}/metrics.csv)\n\n' + table +
                   '\n## Manuscript in progress\n\n'
                   'A preprint-style working draft with sections we will develop together. **No bioRxiv submission or DOI is claimed.**\n\n'
                   '| Section | Working file | Status |\n| --- | --- | --- |\n'
                   '| Abstract | [Abstract](publication/2026-10-06_manuscript/01-abstract.md) | Outline for joint drafting |\n'
                   '| Introduction | [Introduction](publication/2026-10-06_manuscript/02-introduction.md) | Outline for joint drafting |\n'
                   '| Methods | [Methods](publication/2026-10-06_manuscript/03-methods.md) | Existing computational subsection adapted for GitHub |\n'
                   '| Results | [Results](publication/2026-10-06_manuscript/04-results.md) | Source links and writing outline |\n'
                   '| Discussion | [Discussion](publication/2026-10-06_manuscript/05-discussion.md) | Outline and limitations |\n'
                   '| References | [Reference library](publication/2026-10-06_reference-library.md) | PDF locations and retrieval gaps |\n\n'
                   'For scientific work, start with `python3 scripts/workspace.py show`, [campaign overview](overview.md), '
                   '[workspace instructions](docs/WORKSPACE.md), and [challenge requirements](docs/CHALLENGE.md). '
                   'Git includes this selected presentation package; bulk scientific evidence and private bindings remain outside Git.\n')
    (repo / 'README.md').write_text(root_readme)
    print(f'Prepared {len(records)} candidates and {len(panel_rows)} verified species/model pairs under {PANEL}')

if __name__ == '__main__':
    main()
