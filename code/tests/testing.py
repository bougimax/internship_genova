from pathlib import Path
import subprocess
from datetime import datetime
import shutil
import time
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(filename="test.log", level=logging.INFO)


EXECUTABLE_FP = Path(__file__).resolve().parents[1] / "build" / "cdt"
MESHES_FP = Path(__file__).resolve().parents[1] / "Input_file" / "dataset"

meshes = sorted(
    [f for f in MESHES_FP.iterdir() if f.is_file() and f.suffix == ".off"],
    key=lambda f: f.stat().st_size,
)

start_mesh = 1571

num_meshes = len(meshes)

for num_treated_meshes, f in enumerate(meshes[start_mesh:]):
    pre_message = f"Going to treat {MESHES_FP / f.name}, it's {num_treated_meshes+1+start_mesh}/{num_meshes}"
    logger.info(pre_message)
    print(pre_message)
    start = time.time()
    result = subprocess.run(
        [str(EXECUTABLE_FP), "-olvc", f"{MESHES_FP / f.name}"], capture_output=True
    )
    end = time.time()
    if result.stderr != b"":
        error_message = f"Failed on {MESHES_FP / f.name}"
        logger.error(error_message)
        print(error_message)
    succeed_message = f"Treated {MESHES_FP / f.name}, it's {num_treated_meshes+1+start_mesh}/{num_meshes}, it took {end - start:.6f} seconds"
    logger.info(succeed_message)
    print(succeed_message)
