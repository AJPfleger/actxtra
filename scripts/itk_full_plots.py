#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import hist


# Open the ROOT file
file_path = (
    "../../gx2f-push/testexports/itk_output_pixel_1kk/itk_trackstates_fitter.root"
)
rf = uproot.open(file_path)

# # Print the list of keys (objects) in the ROOT file
# print("Keys in the ROOT file:")
# print(rf.keys())

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

# Get first entry for preparing the histograms
df = next(tree.iterate(library="ak", how=dict))

# Get min-max values for eta
data_eta = data_leaf = df["eta_ubs"]
eta_min = min(ak.flatten(data_leaf))
eta_max = max(ak.flatten(data_leaf))

# Create histograms
all_hists = []
for leave in all_leaves:
    data_leaf = df[leave]

    # h_min = min(ak.flatten(data_leaf))
    # h_max = max(ak.flatten(data_leaf))
    h_min = -6
    h_max = 6

    h = hist.Hist(
        hist.axis.Regular(bins=n_bins_leaf, start=h_min, stop=h_max, name=leave),
        hist.axis.Regular(bins=n_bins_eta, start=eta_min, stop=eta_max, name="eta"),
        hist.axis.IntCategory(volume_id_all, name="volume"),
        hist.axis.IntCategory(layer_id_all, name="layer"),
    )

    all_hists.append(h)

# Fill histograms
for df in tree.iterate(library="ak", how=dict):

    # print(df)
    # print("\n")
    data_eta = df["eta_ubs"]
    volume_id = df["volume_id"]
    layer_id = df["layer_id"]

    for leave_i, leave in enumerate(all_leaves):
        data_leaf = df[leave]

        # # TODO how should we handle the surfaces?
        # v_l_pair = (9, 2)
        # log_volume = volume_id == v_l_pair[0]
        # log_layer = layer_id == v_l_pair[1]
        # log_plot = log_volume * log_layer
        # fig, ax = plt.subplots()
        # all_hists[leave_i].fill(ak.flatten(data_leaf[log_plot]), ak.flatten(data_eta[log_plot]))
        all_hists[leave_i].fill(
            ak.flatten(data_leaf),
            ak.flatten(data_eta),
            ak.flatten(volume_id),
            ak.flatten(layer_id),
        )
        # all_hists[leave_i].plot2d(ax=ax, label=f"{leave}", ls="solid")
        # plt.show()

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
#         data_leaf,
#         data_eta,
#         volume_id,
#         layer_id,
#     )

print("*** after filling ***")

volume_layer_pairs = [(9, 2), (9, 4), (16, 2), (16, 4), (16, 6)]
for v_l_pair in volume_layer_pairs:
    v_id = v_l_pair[0]
    l_id = v_l_pair[1]

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
    for leave, ax, h in zip(all_leaves, axes.flatten(), all_hists):
        data_leaf = df[leave]

        h[:, :, v_id * 1j, l_id * 1j].plot2d(ax=ax, label=f"{leave}", ls="solid")

        density = h[:, :, v_id * 1j, l_id * 1j].density()

        bins = h[:, :, v_id * 1j, l_id * 1j].axes[0].edges
        bins = bins[0:-1]
        bins += (bins[1] - bins[0]) / 2

        means_all = []
        rms_all = []
        for density_row in density.T:
            m = sum(density_row * bins) / sum(density_row)
            rms = np.sqrt(sum(density_row * (bins - m) ** 2) / sum(density_row))

            means_all.append(m)
            rms_all.append(rms)

        v_space = -0.285
        for i, (m, rms) in enumerate(zip(means_all, rms_all)):
            s = f"$\\mu$ = {m:.3f}, $\\sigma$ = {rms:.3f}"
            x_pos = -6
            y_pos = 1.35 + i * v_space
            plt.text(x_pos, y_pos, s, fontsize=8, color="w")

    plt.suptitle(f"Volume {v_id}, Layer {l_id}", fontsize=14)
    plt.savefig(f"pull_itk_vol{v_id}_lay{l_id}.pdf")


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
