import matplotlib.pyplot as plt
import numpy as np


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

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    for folder_path in folder_paths:
        MAX_ENERGY_FP = folder_path + "/max_energy.txt"
        max_energy = parse_txt_file(MAX_ENERGY_FP)
        ax.scatter(list(range(len(max_energy))), max_energy, label=folder_path)

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
        ax.scatter(list(range(len(mean_energy))), mean_energy, label=folder_path)
    ax.set_xlabel("Number of iteration of optimization pass")
    ax.set_ylabel("Average energy")
    ax.legend()
    fig.suptitle(
        "Average energy across all the model\n according to the number of optimization pass"
    )
    fig.tight_layout()
    fig.savefig(output_folder + "/mean_energy_evolution.pdf")
    plt.show()

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    for folder_path in folder_paths:
        MEAN_ENERGY_FP = folder_path + "/mean_energy.txt"
        mean_energy = parse_txt_distribution_file(MEAN_ENERGY_FP)
        ax.scatter(list(range(len(mean_energy))), mean_energy, label=folder_path)
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
    fig, ax = plt.subplots()
    distributions = parse_txt_distribution_file(filepath)

    ax.set_xscale("log")
    ax.set_yscale("log")
    for distrib in distributions[::5]:
        counts, bins = np.histogram(distrib, 1000000)
        ax.hist(counts, bins)

    ax.set_xlabel("Energy")
    ax.set_ylabel("Number of tetrahedras")
    ax.legend()
    fig.suptitle("Evolution of energy distribution")
    fig.tight_layout()
    fig.savefig(output_folder + "/energy_distribution_evolution.pdf")
    plt.show()


FOLDER_PATHS = ["2025-06-03_energy_big_float", "2025-06-03_energy_double"]

OUTPUT_FOLDER = "output/test"

DISTRIBUTION_FILE = "test/112544_energy_distribution.txt"

# plot_energy_evolution(FOLDER_PATHS, OUTPUT_FOLDER)

plot_distribution_evolution(DISTRIBUTION_FILE, OUTPUT_FOLDER)
