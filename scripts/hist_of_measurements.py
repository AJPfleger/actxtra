import matplotlib.pyplot as plt
import numpy as np

# Read data
meas_in = []
meas_out = []

filename = "../../gx2f-push/testexports/tmp999_tests"
with open(filename, "r") as file:
    for line in file:
        if line.startswith("QPDATAINPUTMEASUREMENTSSIZE"):
            parts = line.split()
            meas_in.append(int(parts[1]))
        elif line.startswith("QPDATAOUTMEASUREMENTSSIZE"):
            parts = line.split()
            meas_out.append(int(parts[1]))

# Plot the histogram
x = [meas_in, meas_out]
bins = range(np.min(np.concatenate(x)), np.max(np.concatenate(x)) + 2)

left_of_first_bin = np.min(np.concatenate(x)) - 1.0 / 2
right_of_last_bin = np.max(np.concatenate(x)) + 1.0 / 2
bins = np.arange(left_of_first_bin, right_of_last_bin + 1, 1)


fig, ax = plt.subplots()
ax.hist(x, bins=bins, histtype="bar", label=["before fit", "after fit"])
ax.set_title("Histogram of Measurements in the GX2F")
ax.set_xlabel("Number of measurements")
ax.legend(prop={"size": 10})


plt.show()
