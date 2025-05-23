import matplotlib.pyplot as plt


def parse_txt_file(filepath: str) -> list[float]:
    with open(filepath, "r") as file:
        return list(map(float, file.readlines()))


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


FOLDER_PATHS = ["sphere_delta", "sphere_pre_energy"]

OUTPUT_FOLDER = "output/comparison_sphere"

plot_energy_evolution(FOLDER_PATHS, OUTPUT_FOLDER)
