# Alexander J. Pfleger
# 2022-11-25
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case
# params are [y, k]

import numpy as np
import matplotlib.pyplot as plt
import logging

import chi2_utilities as c2u
import propagators


def ai_bi(ri, Vi, xi, params):

    ai = 1 / Vi * np.array([[1, xi], [xi, xi ** 2],])
    bi = ri / Vi * np.array([1, xi])

    return ai, bi


def get_pulls(plot_all, layers=12, cov=0.1):

    ## Initialising
    n_update = 15
    true_params = [12.345, 1]
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    start_params = np.array([0.0, 0.0])  # [y, k]

    # Geometry
    # detector_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5])  # ,10,50,98,99,100])
    detector_layers = np.random.uniform(1, 12, layers)

    # Generate hits
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detector_layers,
        true_params,
        propagators.straight_line_propagator_2D_yk,
        cov,
        True,
    )

    # Fit
    a, updated_params, chi2sum, updated_cov = c2u.gx2f(
        start_params,
        detector_layers,
        cov_meas,
        measurments,
        propagators.straight_line_propagator_2D_yk,
        n_update,
        ai_bi,
    )

    params_res, params_pulls = c2u.calc_res_pulls(
        updated_params, true_params, updated_cov
    )

    y_res, k_res = params_res
    y_pull, k_pull = params_pulls

    if plot_all:
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updated_cov:\n{updated_cov}")
        print(f"pulls: {params_pulls}")
        print("\n")

        max_horizontal = max(detector_layers) + 1
        min_v = min(min(measurments), start_params[0])
        max_v = max(max(measurments), start_params[0])
        delta_v = max_v - min_v
        min_vertical = min_v - 0.3 * delta_v
        max_vertical = max_v + 0.3 * delta_v

        fig, ax = plt.subplots()

        # Detector
        for d in range(len(detector_layers)):
            ax.plot(
                [detector_layers[d], detector_layers[d]],
                [min_vertical, max_vertical],
                "g-",
            )
            ax.plot(detector_layers[d], measurments[d], "gx")

        # Trajectories
        c2u.add_traj_to_plot(
            ax,
            start_params,
            max_horizontal,
            propagators.straight_line_propagator_2D_yk,
            "r",
            "Start Trajectory",
            "-",
        )
        c2u.add_traj_to_plot(
            ax,
            updated_params,
            max_horizontal,
            propagators.straight_line_propagator_2D_yk,
            "b",
            "Final Trajectory",
            "-",
        )
        c2u.add_traj_to_plot(
            ax,
            true_params,
            max_horizontal,
            propagators.straight_line_propagator_2D_yk,
            "k",
            "Unsmeared True Trajectory",
            "-.",
        )

        ax.set(xlabel="x", ylabel="y", title="2D-Fit [y,k]")
        ax.legend()

        fig.savefig("yk-toydetector.pdf")
        plt.show()

    ## root fit
    params_res_root, params_pulls_root = c2u.root_fit(
        detector_layers, measurments, cov, "[0] + [1]*x", true_params
    )
    y_res_root, k_res_root = params_res_root
    y_pull_root, k_pull_root = params_pulls_root

    return y_pull, k_pull, y_res, k_res, y_pull_root, k_pull_root, chi2sum


logging.getLogger().setLevel(logging.INFO)
np.random.seed(10)
draws = 100
layers = 30
bins = int(np.sqrt(draws))
y_pul = []
k_pul = []
y_res_vec = []
k_res_vec = []
y_pul_root = []
k_pul_root = []
chi2sum = []
for d in range(draws):
    print("")  # set this line when using spyder, to make root work correctly
    y_p, k_p, y_res, k_res, y_p_root, k_p_root, c2s = get_pulls(d < 5, layers)
    y_pul.append(y_p)
    k_pul.append(k_p)
    y_res_vec.append(y_res)
    k_res_vec.append(k_res)
    y_pul_root.append(y_p_root)
    k_pul_root.append(k_p_root)
    chi2sum.append(c2s)


c2u.plot_pull_distribution(y_pul, f"y_pulls ({layers} hits)")
c2u.plot_pull_distribution(y_pul_root, f"y_pulls_root ({layers} hits)")
c2u.plot_pull_distribution(k_pul, f"k_pulls ({layers} hits)")
c2u.plot_pull_distribution(k_pul_root, f"k_pulls_root ({layers} hits)")
c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y,k], {layers} hits)")


from matplotlib.patches import Ellipse

cov = np.cov(y_res_vec, k_res_vec)
lambda_, v = np.linalg.eig(cov)
lambda_ = np.sqrt(lambda_)

fig, ax = plt.subplots()
n_scatter_points = 1000
ax.scatter(y_res_vec[:n_scatter_points], k_res_vec[:n_scatter_points])

for j in range(1, 4):
    ell = Ellipse(
        xy=(np.mean(y_res_vec), np.mean(k_res_vec)),
        width=lambda_[0] * j * 2,
        height=lambda_[1] * j * 2,
        angle=np.rad2deg(np.arctan2(*v[:, 0][::-1])),
    )
    ell.set_facecolor("none")
    ell.set_edgecolor("red")
    ax.add_artist(ell)

ax.set(xlabel="y_res", ylabel="k_res", title="y_res vs. k_res, 1e3 runs")
# fig.savefig("yk-y_res-vs-k_res.pdf")
plt.show()
