#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import hist


# Open the ROOT file
file_path = "../../gx2f-push/testexports/itk_output/itk_trackstates_fitter.root"
rf = uproot.open(file_path)

# # Print the list of keys (objects) in the ROOT file
# print("Keys in the ROOT file:")
# print(rf.keys())

# Choose a specific key (object) to extract data from
treename = rf.keys()[-1]  # "trackstates;1"  # Replace with the actual key name
tree = rf[treename]

for df in tree.iterate(library="ak", how=dict):
    print(df["eLOC0_prt"])

    g_x = df["g_x_hit"]
    g_y = df["g_y_hit"]
    g_z = df["g_z_hit"]
    layer_id = df["layer_id"]
    volume_id = df["volume_id"]

    print(np.unique(ak.flatten(layer_id)))
    print(np.unique(ak.flatten(volume_id)))

    g_r = (g_x**2 + g_y**2) ** 0.5

    # plot detector
    log_pix_bar_in = volume_id == 9
    log_pix_bar_out = volume_id == 16
    log_strip_bar = volume_id == 23
    log_plot = log_pix_bar_in + log_pix_bar_out + log_strip_bar

    g_r = g_r[log_plot]
    g_z = g_z[log_plot]
    layer_id = layer_id[log_plot]
    volume_id = volume_id[log_plot]

    fig, ax = plt.subplots()
    scatter = plt.scatter(
        ak.flatten(g_z),
        ak.flatten(g_r),
        # c=ak.flatten(layer_id),
        c=ak.flatten(volume_id),
        cmap=plt.cm.Spectral,
        # vmin=18,
        # vmax=24,
    )
    legend1 = ax.legend(*scatter.legend_elements(), loc="lower left", title="Classes")
    ax.add_artist(legend1)
    plt.show()

    ## Plot pulls and residuals for each layer
    # Create a 2x3 grid of subplots
    fig, axes = plt.subplots(nrows=3, ncols=6, figsize=(12, 8))
    leaves = (
        "res_eLOC0_prt",
        "res_eLOC1_prt",
        "res_ePHI_prt",
        "res_eTHETA_prt",
        "res_eQOP_prt",
        "res_eT_prt",
        "err_eLOC0_prt",
        "err_eLOC1_prt",
        "err_ePHI_prt",
        "err_eTHETA_prt",
        "err_eQOP_prt",
        "err_eT_prt",
        "pull_eLOC0_prt",
        "pull_eLOC1_prt",
        "pull_ePHI_prt",
        "pull_eTHETA_prt",
        "pull_eQOP_prt",
        "pull_eT_prt",
    )

    for i, ax in enumerate(axes.flatten()):
        data_leaf = df[leaves[i]]
        volume_layer_pairs = [(9, 2), (9, 4), (16, 2), (16, 4), (16, 6)]
        # volume_layer_pairs = [(9,2)]

        hmin = min(ak.flatten(data_leaf))
        hmax = max(ak.flatten(data_leaf))
        # print(f"{leaves[i]}: {hmin}, {hmax}")

        # To prevent crash in case of bad distribution
        if hmin == hmax:
            ax.plot([hmin, hmax], [0, 1])
            ax.set(xlabel=leaves[i])
            continue

        for v_l_pair in volume_layer_pairs:
            log_volume = volume_id == v_l_pair[0]
            log_layer = layer_id == v_l_pair[1]
            log_plot = log_volume + log_layer

            # Create a histogram
            h = hist.Hist(
                hist.axis.Regular(bins=100, start=hmin, stop=hmax, name=leaves[i])
            )

            h.fill(ak.flatten(data_leaf[log_plot]))

            h.plot1d(ax=ax, label=f"Vol {v_l_pair[0]}, Lay {v_l_pair[1]}", ls="solid")
            # ax.title.set_text(leaves[i])

    plt.legend()
    plt.show()

    break


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
