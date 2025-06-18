from pathlib import Path
import subprocess
from datetime import datetime
import shutil


EXECUTABLE_FP = Path(__file__).resolve().parents[1] / "build" / "cdt"
MESHES_FP = Path("good_meshes")

ENERGYS_EVOLUTION_FP = Path(
    f"test_{MESHES_FP.name}_{datetime.now().strftime("%Y-%m-%d_%H-%M")}"
)
ENERGYS_EVOLUTION_FP.mkdir()

meshes = sorted(
    [f for f in MESHES_FP.iterdir() if f.is_file() and f.suffix == ".off"],
    key=lambda f: f.stat().st_size,
)

num_meshes = len(meshes)

for num_treated_meshes, f in enumerate(meshes):
    subprocess.run(
        [str(EXECUTABLE_FP), "-olvc", f"{Path.cwd()}/{MESHES_FP.name}/{f.name}"]
    )
    print(f"Treated {f.name}, it's {num_treated_meshes+1}/{num_meshes}")
