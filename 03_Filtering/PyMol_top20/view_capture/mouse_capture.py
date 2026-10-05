from pathlib import Path
from pymol import cmd
VIEW_PATH = str(Path(__file__).with_name('mouse_view.json'))
def save_view():
    import json
    with open(VIEW_PATH, "w") as h: json.dump(list(cmd.get_view()), h)
    print("Saved fixed view to", VIEW_PATH)
    cmd.quit()
cmd.set_key("F5", save_view)
print("Mouse EGFR capture session: set the desired perspective, then press F5.")
