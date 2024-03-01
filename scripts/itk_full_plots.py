#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import hist


def get_bins(edges):
    """
    Calculates the center of the bins from a vector containing edges.
    Such a vector could look like
    h.axes[0].edges
    """
    bins = edges
    bins = bins[0:-1]
    bins += (bins[1] - bins[0]) / 2

    return bins


# Open the ROOT file
file_path = (
    "../../gx2f-push/testexports/itk_output_pixel_1kk/itk_trackstates_fitter.root"
)
rf = uproot.open(file_path)

# Choose a specific key (object) to extract data from
treename = rf.keys()[-1]  # "trackstates;1"  # Replace with the actual key name
tree = rf[treename]

# Surfaces to investigate
volume_layer_pairs = []
volume_layer_pairs.extend([(9, 2), (9, 4), (16, 2), (16, 4), (16, 6)])  # Pixel Barrel
volume_layer_pairs.extend([(23, 2), (23, 4), (23, 6), (23, 8)])  # Strip Barrel
volume_id_all = list(set(v_l_pair[0] for v_l_pair in volume_layer_pairs))
layer_id_all = list(set(v_l_pair[1] for v_l_pair in volume_layer_pairs))

# Binning for the histograms
n_bins_leaf = 100
n_bins_eta = 11

# Interesting leaves
all_leaves = (
    # "res_eLOC0_ubs",
    # "res_eLOC1_ubs",
    # "res_ePHI_ubs",
    # "res_eTHETA_ubs",
    # "res_eQOP_ubs",
    # "res_eT_ubs",
    # "err_eLOC0_ubs",
    # "err_eLOC1_ubs",
    # "err_ePHI_ubs",
    # "err_eTHETA_ubs",
    # "err_eQOP_ubs",
    # "err_eT_ubs",
    "pull_eLOC0_ubs",
    "pull_eLOC1_ubs",
    "pull_ePHI_ubs",
    "pull_eTHETA_ubs",
    "pull_eQOP_ubs",
    # "pull_eT_ubs",
)

# Get min-max values for eta and leaves
eta_max = 1.7
eta_min = -eta_max
h_max = 6
h_min = -h_max

# Create histograms
all_hists = []
for leave in all_leaves:
    h = hist.Hist(
        hist.axis.Regular(bins=n_bins_eta, start=eta_min, stop=eta_max, name="eta"),
        hist.axis.Regular(bins=n_bins_leaf, start=h_min, stop=h_max, name=leave),
        hist.axis.IntCategory(volume_id_all, name="volume"),
        hist.axis.IntCategory(layer_id_all, name="layer"),
    )
    all_hists.append(h)

# Fill histograms
for df in tree.iterate(library="ak", how=dict):
    data_eta = df["eta_ubs"]
    volume_id = df["volume_id"]
    layer_id = df["layer_id"]

    for leave_i, leave in enumerate(all_leaves):
        data_leaf = df[leave]

        all_hists[leave_i].fill(
            ak.flatten(data_eta),
            ak.flatten(data_leaf),
            ak.flatten(volume_id),
            ak.flatten(layer_id),
        )

# # fill with random for faster test
# rand_size = 1_000_000
# data_eta = np.random.normal(size=rand_size) * 2
# volume_id = np.random.randint(20, size=rand_size)
# layer_id = np.random.randint(6, size=rand_size)
#
# for leave_i, leave in enumerate(all_leaves):
#     data_leaf = np.random.normal(size=rand_size) + np.random.normal(size=1)
#
#     all_hists[leave_i].fill(
#         data_eta,
#         data_leaf,
#         volume_id,
#         layer_id,
#     )

print("*** after filling ***")

for v_l_pair in volume_layer_pairs:
    v_id = v_l_pair[0]
    l_id = v_l_pair[1]

    plot_cols = 5
    fig, axes = plt.subplots(nrows=3, ncols=plot_cols, figsize=(12, 8))
    means_all_leaves = []
    rms_all_leaves = []
    for i, (leave, ax, h) in enumerate(
        zip(all_leaves * 3, axes.flatten(), all_hists * 3)
    ):
        if i < plot_cols:
            h[:, :, v_id * 1j, l_id * 1j].plot2d(
                ax=ax, label=f"{leave}", ls="solid", cbar=False
            )

            density = h[:, :, v_id * 1j, l_id * 1j].density()

            bins_eta = get_bins(h[:, :, v_id * 1j, l_id * 1j].axes[0].edges)
            bins_leave = get_bins(h[:, :, v_id * 1j, l_id * 1j].axes[1].edges)

            # Get means for other plots
            means_all_bins = []
            rms_all_bins = []
            for density_row in density:
                m = sum(density_row * bins_leave) / sum(density_row)
                rms = np.sqrt(
                    sum(density_row * (bins_leave - m) ** 2) / sum(density_row)
                )

                means_all_bins.append(m)
                rms_all_bins.append(rms)

            means_all_leaves.append(means_all_bins)
            rms_all_leaves.append(rms_all_bins)

            # Modify plot
            ax.set_title(f"{leave}")
            ax.set_xticks([])
            ax.set_xlabel("")
            if i % plot_cols == 0:
                ax.set_ylabel("pull")
            else:
                ax.set_ylabel("")
                ax.set_yticks([])

        elif i < 2 * plot_cols:
            ax.plot(bins_eta, means_all_leaves[i - plot_cols], ".")
            ax.plot([eta_min, eta_max], [0, 0], "k--")
            ax.set_xlim(eta_min, eta_max)
            ax.set_ylim(-0.1, 0.1)
            ax.set_xticks([])
            ax.set_xlabel("")
            if i % plot_cols == 0:
                ax.set_ylabel("$\\mu$")
            else:
                ax.set_ylabel("")
                ax.set_yticks([])
        else:
            ax.plot(bins_eta, rms_all_leaves[i - 2 * plot_cols], ".")
            ax.plot([eta_min, eta_max], [1, 1], "k--")
            ax.set_xlim(eta_min, eta_max)
            ax.set_ylim(0.8, 1.2)
            ax.set_xlabel("$\\eta$")
            if i % plot_cols == 0:
                ax.set_ylabel("$\\sigma$")
            else:
                ax.set_ylabel("")
                ax.set_yticks([])

    plt.suptitle(f"Volume {v_id}, Layer {l_id}", fontsize=14)
    plt.tight_layout()
    # Adjust spacing between subplots
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.show()
    break
    # plt.savefig(f"pull_itk_vol{v_id}_lay{l_id}.pdf")


# ITk Geometry
# Strip Barrel  Volume 23   Layer 8
#                           Layer 6
#                           Layer 4
#                           Layer 2
#
# Pixel Barrel  Volume 16   Layer 6
#                           Layer 4
#                           Layer 2
#
#               Volume  9   Layer 4
#                           Layer 2
