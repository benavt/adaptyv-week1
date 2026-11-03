"""Render the published native models without prediction or rescoring.

Execute using a Python environment with PyMOL, from the task repository root.
"""
import hashlib
import json
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])
panel = Path('candidates/2026-10-06_readme-panel_001')
manifest = json.loads((panel / 'manifest.json').read_text())
reference = panel / manifest['records'][0]['candidate_id'] / 'Human.cif'
counts = []
for record in manifest['records']:
    folder = panel / record['candidate_id']
    for species in ('Human', 'Mouse'):
        cmd.reinitialize()
        cmd.load(str(reference), 'reference')
        cmd.load(str(folder / (species + '.cif')), 'complex')
        native_count = cmd.count_atoms('complex')
        cmd.align('complex and chain B and name CA', 'reference and chain B and name CA', cycles=0)
        cmd.orient('reference and chain B')
        cmd.rotate('y', -15)
        cmd.rotate('x', 15)
        cmd.delete('reference')
        cmd.hide('everything', 'all')
        cmd.show('cartoon', 'complex and polymer.protein')
        cmd.color('gray80', 'complex and chain B')
        cmd.set_color('binder_cyan', [0.05, 0.66, 0.75])
        cmd.color('binder_cyan', 'complex and chain A')
        if record['score_design_id'] == 'D46_130_fusion':
            cmd.color('magenta', 'complex and chain A and resi 1-14')
        cmd.select('target_histidines', 'complex and chain B and resi 37+100 and resn HIS')
        assert cmd.count_atoms('target_histidines and name CA') == 2
        cmd.show('sticks', 'target_histidines')
        cmd.color('forest', 'target_histidines')
        cmd.show('sticks', 'complex and not polymer.protein')
        cmd.set('cartoon_fancy_helices', 1)
        cmd.set('cartoon_loop_radius', 0.22)
        cmd.set('stick_radius', 0.18)
        cmd.set('orthoscopic', 1)
        cmd.set('antialias', 2)
        cmd.set('ray_shadows', 0)
        cmd.set('ambient', 0.55)
        cmd.set('specular', 0.15)
        cmd.set('ray_opaque_background', 1)
        cmd.bg_color('white')
        cmd.viewport(720, 520)
        cmd.zoom('complex', buffer=1, complete=0)
        cmd.deselect()
        cmd.png(str(folder / (species + '.png')), width=720, height=520, dpi=150, ray=1)
        cmd.save(str(folder / (species + '.pse')))
        cmd.reinitialize()
        cmd.load(str(folder / (species + '.pse')))
        assert cmd.count_atoms('complex') == native_count
        counts.append({'candidate_id': record['candidate_id'], 'species': species,
                       'native_and_reopened_session_atom_count': native_count})
        print(record['candidate_id'], species, flush=True)
manifest['rendering'] = {'created_at': datetime.now(ZoneInfo('America/New_York')).isoformat(),
                         'engine': 'PyMOL', 'version': list(cmd.get_version()),
                         'reference_model': str(reference),
                         'camera': 'EGFR CA alignment to first panel human model, no outlier rejection; reference receptor oriented, y -15 and x +15 degrees; whole complex framed',
                         'coordinate_policy': 'Native CIFs unchanged. PSE coordinates rigidly transformed for display; no score recalculation.',
                         'PNG_size': [720, 520], 'ray': True, 'reopened_session_counts': counts}
manifest['outputs'] = {str(p.relative_to(panel)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in panel.rglob('*') if p.is_file() and p.name != 'manifest.json'}
(panel / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
cmd.quit()
