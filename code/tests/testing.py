from pathlib import Path
import subprocess
from datetime import datetime
import shutil
import time
import logging
import os
import numpy as np
from math import inf


def parse_node_file(path):
    nodes = []
    with open(path, "r") as f:
        for line in f.readlines()[1:-1]:
            nodes.append(np.array(list(map(float, line.split()[1:]))))
    return nodes


def parse_ele_file(path):
    tets = []
    with open(path, "r") as f:
        for line in f.readlines()[1:-1]:
            tets.append(tuple(map(int, line.split()[1:])))
    return tets


def compute_energy_tetrahedra(p1, p2, p3, p4):
    j1 = 0.5 * (p1 + p3 - p2 - p4)
    j2 = 0.5 * (p1 - p3 - p2 + p4)
    j3 = 0.5 * (p3 - p1 - p2 + p4)

    J = np.array([j1, j2, j3])

    return np.trace(np.matmul(J, np.transpose(J))) / ((np.linalg.det(J)) ** (2 / 3))


def compute_energys(model_name):

    nodes = parse_node_file(f"{MESHES_FP}/{model_name}/{model_name}.1.node")
    tets = parse_ele_file(f"{MESHES_FP}/{model_name}/{model_name}.1.ele")

    max_energy, mean_energy = 0, 0
    n_tets = len(tets)

    for i1, i2, i3, i4 in tets:
        e = compute_energy_tetrahedra(nodes[i1], nodes[i2], nodes[i3], nodes[i4])
        max_energy = max(max_energy, e)
        mean_energy += e / n_tets

    return (round(max_energy, 5), round(mean_energy, 5))


logger = logging.getLogger(__name__)
logging.basicConfig(filename="test.log", level=logging.INFO)


EXECUTABLE_FP = Path(__file__).resolve().parents[1] / "build" / "cdt"
MESHES_FP = Path(__file__).resolve().parents[1] / "Input_file" / "dataset"
TETGEN_FP = Path(__file__).resolve().parents[0] / "TetGen" / "tetgen"

meshes = sorted(
    [f for f in MESHES_FP.iterdir() if f.is_file() and f.suffix == ".off"],
    key=lambda f: f.stat().st_size,
)

start_mesh = 0

num_meshes = len(meshes)

for num_treated_meshes, f in enumerate(meshes[start_mesh:]):

    # OUR PART
    model_name = f.stem

    pre_message = f"Going to treat {model_name}, it's {num_treated_meshes+1+start_mesh}/{num_meshes}"
    logger.info(pre_message)
    print(pre_message)
    start = time.time()
    result = subprocess.run(
        [str(EXECUTABLE_FP), "-olvc", f"{MESHES_FP / f.name}"], capture_output=True
    )
    end = time.time()
    if result.stderr != b"":
        error_message = f"Failed on {model_name}"
        logger.error(error_message)
        print(error_message)
    else:
        succeed_message = f"Treated {model_name}, it's {num_treated_meshes+1+start_mesh}/{num_meshes}, it took {end - start:.6f} seconds"
        logger.info(succeed_message)
        print(succeed_message)

    # TETGEN PART

    pre_message = f"Launching tetgen on {model_name}"
    logger.info(pre_message)
    print(pre_message)
    start = time.time()
    try:
        result_tetgen = subprocess.run(
            [str(TETGEN_FP), "-O2/7/3", "-pqFBQ", f"{MESHES_FP / f.name}"],
            timeout=300,
            capture_output=True,
        )
        end = time.time()
        if result_tetgen.returncode == 0:
            succeed_message = (
                f"Tetgen treated {model_name}, it took {end - start:.6f} seconds"
            )
            logger.info(succeed_message)
            print(succeed_message)
            os.remove(f"{MESHES_FP}/{model_name}.1.smesh")
            for ext in ["ele", "node"]:
                shutil.move(
                    f"{MESHES_FP}/{model_name}.1.{ext}",
                    f"{MESHES_FP}/{model_name}/{model_name}.1.{ext}",
                )
            with open(Path(f"{MESHES_FP}/{model_name}/result_tetgen.txt"), "w") as f:
                max_energy, mean_energy = compute_energys(model_name)
                f.write(f"max_energy : {round(float(max_energy),5)}\n")
                f.write(f"mean_energy : {round(float(mean_energy),5)}\n")
                f.write(f"time : {end - start:.6f}s")
                with open(
                    Path(f"{MESHES_FP}/{model_name}/mean_energy.txt"), "r"
                ) as f_cdt:
                    mean_energy_cdt = float(f_cdt.readlines()[-1].strip())
                print(
                    f"CDT mean energy is {mean_energy_cdt}, while TetGen's is {mean_energy}, best is {"CDT" if mean_energy_cdt < mean_energy else ("TetGen" if mean_energy_cdt > mean_energy else "both")}"
                )
                with open(
                    Path(f"{MESHES_FP}/{model_name}/max_energy.txt"), "r"
                ) as f_cdt:
                    max_energy_cdt = float(f_cdt.readlines()[-1].strip())
                print(
                    f"CDT max energy is {max_energy_cdt}, while TetGen's is {max_energy}, best is {"CDT" if max_energy_cdt < max_energy else ("TetGen" if max_energy_cdt > max_energy else "both")}"
                )

        else:
            error_message = f"Tetgen failed on {model_name}"
            logger.error(error_message)
            print(error_message)
            with open(Path(f"{MESHES_FP}/{model_name}/result_tetgen.txt"), "w") as f:
                max_energy, mean_energy = compute_energys(model_name)
                f.write(f"max_energy : {inf}\n")
                f.write(f"mean_energy : {inf}\n")
                f.write(f"time : {inf}")

    except subprocess.TimeoutExpired:
        error_message = f"Tetgen failed on {model_name}"
        logger.error(error_message)
        print(error_message)
        with open(Path(f"{MESHES_FP}/{model_name}/result_tetgen.txt"), "w") as f:
            f.write(f"max_energy : {inf}\n")
            f.write(f"mean_energy : {inf}\n")
            f.write(f"time : {inf}")
