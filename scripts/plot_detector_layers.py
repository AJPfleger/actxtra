#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np


def plot_detector(df):
    g_x = df["g_x_hit"]
    g_y = df["g_y_hit"]
    g_z = df["g_z_hit"]
    volume_id = df["volume_id"]

    # print(np.unique(ak.flatten(layer_id)))
    # print(np.unique(ak.flatten(volume_id)))

    g_r = (g_x**2 + g_y**2) ** 0.5

    fig, ax = plt.subplots()
    scatter = plt.scatter(
        ak.flatten(g_z),
        ak.flatten(g_r),
        # c=ak.flatten(layer_id),
        c=ak.flatten(volume_id),
        cmap="Spectral",
        s=1
    )
    # title
    plt.title("ODD Detector " + file_path)

    regions = []
    # print("ids", np.unique(ak.flatten(volume_id)))
    for i_v_id in np.unique(ak.flatten(volume_id)):
        # get x and y max and min value, filter out nan with ak
        x_min = ak.min(ak.flatten(g_z[volume_id == i_v_id][~ak.is_none(g_z[volume_id == i_v_id])]))
        x_max = ak.max(ak.flatten(g_z[volume_id == i_v_id][~ak.is_none(g_z[volume_id == i_v_id])]))
        y_min = ak.min(ak.flatten(g_r[volume_id == i_v_id][~ak.is_none(g_r[volume_id == i_v_id])]))
        y_max = ak.max(ak.flatten(g_r[volume_id == i_v_id][~ak.is_none(g_r[volume_id == i_v_id])]))

        # if x or y values are infinity
        if np.isinf(x_min) or np.isinf(x_max) or np.isinf(y_min) or np.isinf(y_max):
            continue
        x_mid = ((x_max - x_min) / 3) + x_min
        y_mid = ((y_max - y_min) / 2) + y_min

        # print("vol", i_v_id, x_min, x_max, y_min, y_max, x_mid, y_mid)

        regions.append({"x": x_mid, "y": y_mid, "label": f"vol {i_v_id}"})


    for region in regions:
        ax.annotate(
            region["label"],
            xy=(region["x"], region["y"]),
            xytext=(region["x"], region["y"])
        )
    plt.show()


# Open the ROOT file
file_path = "trackstates_3eta_1kkevents.root"
rf = uproot.open(file_path)

# Choose a specific key (object) to extract data from
treename = rf.keys()[-1]  # "trackstates;1"  # Replace with the actual key name
tree = rf[treename]

for df in tree.iterate(library="ak", how=dict):
    plot_detector(df)

    break
