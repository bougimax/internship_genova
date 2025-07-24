import matplotlib.pyplot as plt
import numpy as np
import imageio
from pathlib import Path
import os


def parse_txt_file(filepath: str) -> list[float]:
    with open(filepath, "r") as file:
        return list(map(float, file.readlines()))


def parse_txt_distribution_file(filepath: str) -> list[list[float]]:
    output = []
    with open(filepath, "r") as file:
        for line in file.readlines():
            output.append(list(map(float, line.split())))
    return output


def plot_energy_evolution(folder_paths: list[str], output_folder: str):
    if not Path(output_folder).exists():
        os.mkdir(output_folder)

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    for folder_path in folder_paths:
        MAX_ENERGY_FP = folder_path + "/max_energy.txt"
        max_energy = parse_txt_file(MAX_ENERGY_FP)
        ax.scatter(
            [i / 6 for i in range(len(max_energy))], max_energy, label=folder_path
        )

    ax.set_xlabel("Number of iteration of optimization pass")
    ax.set_ylabel("Maximal energy")
    ax.legend()
    fig.suptitle(
        "Maximal energy across all the model\n according to the number of optimization pass"
    )
    fig.tight_layout()
    fig.savefig(output_folder + "/max_energy_evolution.pdf")

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    for folder_path in folder_paths:
        MEAN_ENERGY_FP = folder_path + "/mean_energy.txt"
        mean_energy = parse_txt_file(MEAN_ENERGY_FP)
        ax.scatter(
            [i / 6 for i in range(len(mean_energy))], mean_energy, label=folder_path
        )
    ax.set_xlabel("Number of iteration of optimization pass")
    ax.set_ylabel("Average energy")
    ax.legend()
    fig.suptitle(
        "Average energy across all the model\n according to the number of optimization pass"
    )
    fig.tight_layout()
    fig.savefig(output_folder + "/mean_energy_evolution.pdf")
    plt.show()


def plot_distribution_evolution(filepath: str, output_folder: str):
    distributions = parse_txt_distribution_file(filepath)

    max_range = 1000

    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlim(3, max_range)
    ax.set_xlabel("Energy")
    ax.set_ylabel("Number of tetrahedras")
    fig.suptitle("Evolution of energy distribution")
    fig.tight_layout()

    for i, distrib in zip(
        ["last", "fifth", "initial"],
        [distributions[-1], distributions[5], distributions[0]],
    ):
        distrib = [x for x in distrib if x <= max_range]
        ax.hist(distrib, bins=5000, label=f"{i} distribution", density=True)

    ax.legend()
    fig.savefig(output_folder + f"/energy_distribution_evolution.pdf")

    for i, distrib in enumerate(distributions):
        distrib = [x for x in distrib if x <= max_range]
        fig, ax = plt.subplots()
        ax.set_xscale("log")
        ax.set_ylim(0, 1)
        ax.set_xlim(3, max_range)
        ax.hist(distrib, bins=5000, label=f"{i}th distribution", density=True)
        ax.legend()
        ax.set_xlabel("Energy")
        ax.set_ylabel("Number of tetrahedras")
        fig.suptitle("Evolution of energy distribution")
        fig.tight_layout()
        fig.savefig(output_folder + f"/energy_distribution_evolution_{i}.png")
        fig.savefig(output_folder + f"/energy_distribution_evolution_{i}.pdf")
        print(f"Done {i}th pass")

    frames = []

    for i in range(len(distributions)):
        image = imageio.v2.imread(
            output_folder + f"/energy_distribution_evolution_{i}.png"
        )
        frames.append(image)

    imageio.mimsave(
        output_folder + "/energy_distribution_evolution.gif",
        frames,
        fps=5,
        loop=1,
    )


FOLDER_PATHS = ["112544"]

OUTPUT_FOLDER = "output/energy_evolution_112544"

DISTRIBUTION_FILE = "112544/energy_distribution.txt"

plot_energy_evolution(FOLDER_PATHS, OUTPUT_FOLDER)

plot_distribution_evolution(DISTRIBUTION_FILE, OUTPUT_FOLDER)
