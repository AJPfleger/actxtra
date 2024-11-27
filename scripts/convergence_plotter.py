import re
import matplotlib.pyplot as plt
import numpy as np


def extract_params_from_log(file_path, n_params):
    with open(file_path, "r") as file:
        log_data = file.readlines()

    pattern = re.compile(r"deltaParamsExtended:")  # To detect start of each section
    results = []
    current_params = []

    for line in log_data:
        if "deltaParamsExtended:" in line:
            if current_params:
                results.append(current_params)  # Save the previous set
            current_params = []  # Start a new set

        elif len(current_params) < n_params:
            try:
                # Try to convert the line to a float
                current_params.append(float(line.strip()))
            except ValueError:
                # If conversion fails, skip this line
                continue

    # Save the last occurrence
    if current_params:
        results.append(current_params)

    return results


def extract_to_array(jagged_list):
    n_params = 0
    for vec in jagged_list:
        n_params = max(n_params, len(vec))
        print(f"size = {len(vec)}\tn_params = {n_params}\n{vec}")

    steps = len(jagged_list)

    dejagged_array = np.empty([n_params, steps]) * np.nan

    for i, vec in enumerate(jagged_list):
        for j, param in enumerate(vec):
            dejagged_array[j, i] = param
        # n_params = max(n_params, len(vec))

    print(dejagged_array[:, 1:])

    # we need to remove the first column, since it is empty
    return dejagged_array[:, 1:]


def plot_params_over_time(occurrences):
    # occurrences should be in the form of [n_steps, n_params]
    occurrences = np.array(occurrences)
    n_params, n_steps = occurrences.shape
    # n_steps = 30
    time_steps = range(n_steps)

    plt.figure(figsize=(10, 8))

    # Top subplot for bound params (first 6 parameters)
    plt.subplot(2, 1, 1)
    for i in range(6):
        plt.plot(time_steps, occurrences[i, time_steps], label=f"Param {i+1}")
    plt.title("Bound Params")
    plt.xlabel("nUpdate")
    plt.ylabel("delta parameter")
    plt.legend(loc="upper right", fontsize="small")
    plt.xlim([0, n_steps - 1])

    if n_params < 6:
        return

    # Bottom subplot for scattering angles (rest of the parameters)
    plt.subplot(2, 1, 2)
    for i in range(6, n_params):
        plt.plot(time_steps, occurrences[i, time_steps], label=f"Param {i+1}")
    plt.title("Scattering Angles")
    plt.xlabel("nUpdate")
    plt.ylabel("delta parameter")
    plt.legend(loc="upper right", fontsize="small")
    plt.xlim([0, n_steps - 1])

    plt.tight_layout()
    plt.show()


# Usage
log_file_path = "<full_path>/physmon/logs/dev_trackfitting_gx2f_vs_kf.log"
n_params = 14  # Set this to the correct number of params
occurrences = extract_params_from_log(log_file_path, n_params)

# Now `occurrences` contains all the collected arrays of deltaParamsExtended
print(occurrences)

dejagged = extract_to_array(occurrences)

# Usage
plot_params_over_time(dejagged)
