import re
import numpy as np
import matplotlib.pyplot as plt

filename = "grid_search_stds.txt"

# Read data from the file
with open(filename, "r") as file:
    lines = file.readlines()

# Regular expression to extract floats from the given data
float_pattern = r"-?\d+\.\d+"

# Extracting ETA, PHI, and x1 values
data = []
for line in lines:
    matches = re.findall(r"-?\d+\.\d+", line)
    print(matches)
    if matches:
        data.append(matches)


# Unpacking the matches
eta_values, phi_values, x1_values, x2_values, x3_values, x4_values, x5_values = zip(
    *data
)

# Converting strings to floats
eta_values = list(map(float, eta_values))
phi_values = list(map(float, phi_values))
phi_values = np.degrees(np.array(phi_values))
x1_values = list(map(float, x1_values))
x2_values = list(map(float, x2_values))
x3_values = list(map(float, x3_values))
x4_values = list(map(float, x4_values))
x5_values = list(map(float, x5_values))


# Create subplots for x1 to x5
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

# Plot scatter for ETA and PHI with x1
sc1 = axs[0, 0].scatter(eta_values, phi_values, c=x1_values, cmap="viridis", marker="o")
axs[0, 0].set_xlabel("$\eta$")
axs[0, 0].set_ylabel("$\phi$")
axs[0, 0].set_title("pull_eLOC0_fit")

# Plot scatter for ETA and PHI with x2
sc2 = axs[0, 1].scatter(eta_values, phi_values, c=x2_values, cmap="viridis", marker="o")
axs[0, 1].set_xlabel("$\eta$")
axs[0, 1].set_ylabel("$\phi$")
axs[0, 1].set_title("pull_eLOC1_fit")

# Plot scatter for ETA and PHI with x3
sc3 = axs[0, 2].scatter(eta_values, phi_values, c=x3_values, cmap="viridis", marker="o")
axs[0, 2].set_xlabel("$\eta$")
axs[0, 2].set_ylabel("$\phi$")
axs[0, 2].set_title("pull_ePHI_fit")

# Plot scatter for ETA and PHI with x4
sc4 = axs[1, 0].scatter(eta_values, phi_values, c=x4_values, cmap="viridis", marker="o")
axs[1, 0].set_xlabel("$\eta$")
axs[1, 0].set_ylabel("$\phi$")
axs[1, 0].set_title("pull_eTHETA_fit")

# Plot scatter for ETA and PHI with x5
sc5 = axs[1, 1].scatter(eta_values, phi_values, c=x5_values, cmap="viridis", marker="o")
axs[1, 1].set_xlabel("$\eta$")
axs[1, 1].set_ylabel("$\phi$")
axs[1, 1].set_title("pull_eQOP_fit")

# Hide the empty subplot
axs[1, 2].axis("off")

# Add colorbars to the right of the subplots
cbar1 = fig.colorbar(sc1, ax=axs[0, 0], pad=0.1)
cbar2 = fig.colorbar(sc2, ax=axs[0, 1], pad=0.1)
cbar3 = fig.colorbar(sc3, ax=axs[0, 2], pad=0.1)
cbar4 = fig.colorbar(sc4, ax=axs[1, 0], pad=0.1)
cbar5 = fig.colorbar(sc5, ax=axs[1, 1], pad=0.1)

plt.tight_layout()
# plt.show()
plt.savefig("2024-01-24_stds.pdf")
