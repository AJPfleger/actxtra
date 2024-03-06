import matplotlib.pyplot as plt
import numpy as np


def flatten2(list_of_list_of_lists):
    step1 = [x for xs in list_of_list_of_lists for x in xs]
    return [x for xs in step1 for x in xs]


energies = [1, 10, 100]

# Read data
meas_io_all = []
for e in energies:
    filename = f"../../gx2f-push/testexports/tmp999_tests{e}"
    meas_in = []
    meas_out = []
    with open(filename, "r") as file:
        for line in file:
            if line.startswith("QPDATAINPUTMEASUREMENTSSIZE"):
                parts = line.split()
                meas_in.append(int(parts[1]))
            elif line.startswith("QPDATAOUTMEASUREMENTSSIZE"):
                parts = line.split()
                meas_out.append(int(parts[1]))
    # meas_out = meas_out * 3
    meas_io_all.append([meas_in, meas_out])

# Prepare histogram binning
left_of_first_bin = min(flatten2(meas_io_all)) - 1.0 / 2
right_of_last_bin = max(flatten2(meas_io_all)) + 1.0 / 2
bins = np.arange(left_of_first_bin, right_of_last_bin + 1, 1)

# Plot histograms
fig, axs = plt.subplots(1, len(energies), figsize=(15, 5))

for i, meas_io in enumerate(meas_io_all):
    ax = axs[i]
    ax.hist(meas_io, bins=bins, histtype="bar", label=["before fit", "after fit"])
    ax.set_title(f"Histogram of Measurements for {energies[i]} GeV")
    ax.set_xlabel("Number of measurements")
    ax.legend(prop={"size": 10})

plt.tight_layout()
# plt.show()
plt.savefig("hist_of_measurements_100k_events.png")
