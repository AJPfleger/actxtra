#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np

# Open the ROOT file
file_path = "../../gx2f-push/testexports/itk_output/itk_trackstates_fitter.root"
rf = uproot.open(file_path)

# # Print the list of keys (objects) in the ROOT file
# print("Keys in the ROOT file:")
# print(rf.keys())

# Choose a specific key (object) to extract data from
treename = rf.keys()[-1] #"trackstates;1"  # Replace with the actual key name
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
    legend1 = ax.legend(*scatter.legend_elements(),
                        loc="lower left", title="Classes")
    ax.add_artist(legend1)
    plt.show()



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
