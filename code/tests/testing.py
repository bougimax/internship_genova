from pathlib import Path
import random as rd
import subprocess
from datetime import datetime
import shutil
import time
import logging
import os
import numpy as np
from math import inf
import statistics


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
logging.basicConfig(
    filename="test.log",
    level=logging.INFO,
    format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


EXECUTABLE_FP = Path(__file__).resolve().parents[1] / "build" / "cdt"
MESHES_FP = Path(__file__).resolve().parents[1] / "Input_file"
TETGEN_FP = Path(__file__).resolve().parents[0] / "TetGen" / "tetgen"

meshes = sorted(
    [f for f in MESHES_FP.iterdir() if f.is_file() and f.suffix == ".off"],
    key=lambda f: f.stat().st_size,
)

start_mesh = 0

num_meshes = len(meshes)

tetgen_better_on_mean, tetgen_better_on_max, tetgen_better_on_time = 0, 0, 0
cdt_better_on_mean, cdt_better_on_max, cdt_better_on_time = 0, 0, 0
tetgen_fails = 0
num_both_treated = 0
relative_error_mean = []
relative_error_max = []
average_relative_error_mean = 0
average_relative_error_max = 0
tie_max, tie_mean, tie_duration = 0, 0, 0

# test_meshes = rd.choices(meshes[:1000], k=100)
test_meshes = meshes[:1000]

for num_treated_meshes, f in enumerate(test_meshes[start_mesh:]):

    # OUR PART
    model_name = f.stem

    pre_message = (
        f"Going to treat {model_name}, {num_treated_meshes+1+start_mesh}/{num_meshes}"
    )
    logger.info(pre_message)
    print(pre_message)
    start = time.time()
    result = subprocess.run(
        [str(EXECUTABLE_FP), "-olv", f"{MESHES_FP / f.name}"], capture_output=True
    )
    end = time.time()
    duration_cdt = end - start
    if result.stderr != b"":
        error_message = f"Failed on {model_name}\n"
        logger.error(error_message)
        print(error_message)
    else:
        succeed_message = f"Treated {model_name}, it took {end - start:.6f} seconds\n"
        logger.info(succeed_message)
        print(succeed_message)

    # TETGEN PART

    pre_message = f"Launching tetgen on {model_name}"
    logger.info(pre_message)
    print(pre_message)
    start = time.time()
    try:
        result_tetgen = subprocess.run(
            [str(TETGEN_FP), "-pqFBQT0", f"{MESHES_FP / f.name}"],
            timeout=10,
            capture_output=True,
        )
        end = time.time()
        duration_tetgen = end - start
        if result_tetgen.returncode == 0:
            breaking = False
            num_both_treated += 1
            succeed_message = (
                f"Tetgen treated {model_name}, it took {end - start:.6f} seconds\n"
            )
            logger.info(succeed_message)
            print(succeed_message)
            os.remove(f"{MESHES_FP}/{model_name}.1.smesh")
            for ext in ["ele", "node"]:
                if not Path(f"{MESHES_FP}/{model_name}.1.{ext}").exists():
                    breaking = True
                    break
                shutil.move(
                    f"{MESHES_FP}/{model_name}.1.{ext}",
                    f"{MESHES_FP}/{model_name}/{model_name}.1.{ext}",
                )
            if breaking:
                tetgen_fails += 1
                error_message = f"Tetgen failed on {model_name}\n"
                logger.error(error_message)
                print(error_message)
                with open(
                    Path(f"{MESHES_FP}/{model_name}/result_tetgen.txt"), "w"
                ) as f:
                    f.write(f"max_energy : {inf}\n")
                    f.write(f"mean_energy : {inf}\n")
                    f.write(f"time : {inf}")
                continue
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
                relative_error_mean.append(
                    (mean_energy_cdt - mean_energy) / mean_energy
                )
                relative_error_max.append((max_energy_cdt - max_energy) / max_energy)
                average_relative_error_mean = (
                    (average_relative_error_mean * (num_both_treated - 1))
                    + ((mean_energy_cdt - mean_energy) / mean_energy)
                ) / num_both_treated
                average_relative_error_max = (
                    (average_relative_error_max * (num_both_treated - 1))
                    + ((max_energy_cdt - max_energy) / max_energy)
                ) / num_both_treated
                tetgen_better_on_max += max_energy < max_energy_cdt
                cdt_better_on_max += max_energy > max_energy_cdt
                tetgen_better_on_mean += mean_energy < mean_energy_cdt
                cdt_better_on_mean += mean_energy > mean_energy_cdt
                cdt_better_on_time += duration_cdt < duration_tetgen
                tetgen_better_on_time += duration_cdt > duration_tetgen
                tie_max += max_energy == max_energy_cdt
                tie_mean += mean_energy == mean_energy_cdt
                tie_duration += duration_cdt == duration_tetgen
                print(
                    f"CDT max energy is {max_energy_cdt}, while TetGen's is {max_energy}, best is {"CDT" if max_energy_cdt < max_energy else ("TetGen" if max_energy_cdt > max_energy else "both")}\n"
                )

        else:
            tetgen_fails += 1
            error_message = f"Tetgen failed on {model_name}\n"
            logger.error(error_message)
            print(error_message)
            with open(Path(f"{MESHES_FP}/{model_name}/result_tetgen.txt"), "w") as f:
                f.write(f"max_energy : {inf}\n")
                f.write(f"mean_energy : {inf}\n")
                f.write(f"time : {inf}")

    except subprocess.TimeoutExpired:
        tetgen_fails += 1
        error_message = f"Tetgen failed on {model_name}\n"
        logger.error(error_message)
        print(error_message)
        with open(Path(f"{MESHES_FP}/{model_name}/result_tetgen.txt"), "w") as f:
            f.write(f"max_energy : {inf}\n")
            f.write(f"mean_energy : {inf}\n")
            f.write(f"time : {inf}")

    print(
        f"""Cdt won: {cdt_better_on_max}/{num_treated_meshes + 1 - tetgen_fails - tie_max} on max energy
         {cdt_better_on_mean}/{num_treated_meshes + 1- tetgen_fails - tie_mean} on mean energy
         {cdt_better_on_time}/{num_treated_meshes + 1- tetgen_fails - tie_duration} on time
Tetgen failed {tetgen_fails}/{num_treated_meshes+1}\n -----------------------------------------------"""
    )
