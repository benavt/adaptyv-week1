from pymol import cmd
from pathlib import Path
import json, os, subprocess

species = os.environ.get('EGFR_SPECIES', '').capitalize()
if species not in ('Mouse', 'Human'):
    raise SystemExit('Set EGFR_SPECIES=Mouse or EGFR_SPECIES=Human')
root = Path.cwd() / '01_Staging' / 'MD' / 'EGFR_D3' / species
frames = sorted((root / 'md_10ns' / 'processed' / 'pdb_frames').glob('frame*.pdb'), key=lambda p: int(p.stem[5:]))
view_path = root / f'egfr_d3_{species.lower()}_view.json'
if not view_path.exists():
    raise SystemExit(f'Missing {view_path}; save the PyMOL view first')
view = json.loads(view_path.read_text())
obj = f'egfr_{species.lower()}'
cmd.reinitialize()
cmd.load(str(frames[0]), obj)
cmd.hide('everything', obj)
cmd.show('cartoon', f'{obj} and polymer.protein')
base_color = 'skyblue' if species == 'Mouse' else 'lightorange'
cmd.color(base_color, f'{obj} and polymer.protein')
cmd.set_color('his_purple', [0.60, 0.00, 0.80])
cmd.select('egfr_all_histidines', f'{obj} and polymer.protein and (resn HIS or resn HISH)')
cmd.color('his_purple', 'egfr_all_histidines')
if species == 'Mouse':
    cmd.color('red', f'{obj} and polymer.protein and resn LYS+ARG')
    cmd.color('blue', f'{obj} and polymer.protein and resn ASP+GLU')
else:
    cmd.color('blue', f'{obj} and polymer.protein and resn LYS+ARG')
    cmd.color('red', f'{obj} and polymer.protein and resn ASP+GLU')
cmd.show('sticks', 'egfr_all_histidines')
cmd.show('sticks', f'{obj} and polymer.protein and resn LYS+ARG+ASP+GLU')
cmd.set('stick_radius', 0.18)
cmd.bg_color('white')
cmd.set('ray_opaque_background', 0)
cmd.set_view(view)
cmd.mset(f'1-{len(frames)}')
for i, frame_path in enumerate(frames, 1):
    cmd.load(str(frame_path), obj, state=i)
    cmd.frame(i)
image_dir = root / 'movie_frames'
image_dir.mkdir(exist_ok=True)
for i in range(1, len(frames) + 1):
    cmd.frame(i)
    cmd.png(str(image_dir / f'frame{i:04d}.png'), width=800, height=600, ray=0, quiet=1)
ffmpeg = '/opt/homebrew/Caskroom/miniforge/base/envs/thermoMPNN/bin/ffmpeg'
output = root / f'egfr_d3_{species.lower()}_10ns.mp4'
env = dict(os.environ, DYLD_LIBRARY_PATH='/opt/homebrew/Caskroom/miniforge/base/envs/thermoMPNN/lib')
subprocess.run([ffmpeg, '-y', '-framerate', '30', '-i', str(image_dir / 'frame%04d.png'), '-c:v', 'mpeg4', '-q:v', '3', '-pix_fmt', 'yuv420p', str(output)], env=env, check=True)
cmd.quit()
