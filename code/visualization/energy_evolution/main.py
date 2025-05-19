import matplotlib.pyplot as plt


def parse_txt_file(filepath: str) -> list[float]:
    with open(filepath, "r") as file:
        return list(map(float, file.readlines()))


def plot_energy_evolution(folder_path: str):
    MAX_ENERGY_FP = folder_path + "/max_energy.txt"
    MEAN_ENERGY_FP = folder_path + "/mean_energy.txt"
    max_energy = parse_txt_file(MAX_ENERGY_FP)
    mean_energy = parse_txt_file(MEAN_ENERGY_FP)

    fig, ax = plt.subplots()
    ax.scatter(list(range(len(max_energy))), max_energy)
    ax.set_xlabel("Number of iteration of optimization pass")
    ax.set_ylabel("Maximal energy")
    fig.suptitle(
        "Maximal energy across all the model\n according to the number of optimization pass"
    )
    fig.tight_layout()
    fig.savefig(folder_path + "/max_energy_evolution.pdf")

    fig, ax = plt.subplots()
    ax.scatter(list(range(len(mean_energy))), mean_energy)
    ax.set_xlabel("Number of iteration of optimization pass")
    ax.set_ylabel("Average energy")
    fig.suptitle(
        "Average energy across all the model\n according to the number of optimization pass"
    )
    fig.tight_layout()
    fig.savefig(folder_path + "/mean_energy_evolution.pdf")


FOLDER_PATH = "2025-05-19"
plot_energy_evolution(FOLDER_PATH)
