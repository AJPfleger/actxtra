# Alexander J. Pfleger
# 2022-01-24
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 1D case with just a single degree of freedom

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

import chi2_utilities as c2u
import propagators


def get_pulls(plot_all, layers=12, cov=0.1):

    ## Initialising
    n_update = 2

    # Parameter
    start_params = 0.0  # [y0]
    delta_params = 0

    # Geometry
    # detector_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5,10,10.1,50,98,99,100])
    detector_layers = np.random.uniform(1, 100, layers)

    # Hits
    true_params = 12.345
    # true_params = [np.random.uniform(-9,9)]
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detector_layers, true_params, propagators.straight_line_propagator_1D, cov
    )

    updated_params = start_params

    ## Iterating and updating parameters
    for _ in range(n_update):

        updated_params = updated_params + delta_params

        # Iterate over surfaces
        a = 0
        b = 0
        chi2sum = 0
        for d in range(len(detector_layers)):
            x = detector_layers[d]
            Vi = cov_meas[d]
            ri = measurments[d] - propagators.straight_line_propagator_1D(
                updated_params, x
            )
            chi2sum += c2u.chi2_1D(Vi, ri)

            ai = 1 / Vi
            bi = ri / Vi

            a += ai
            b += bi

        delta_params = b / a

    updated_cov = 1 / a

    y_res = updated_params - true_params
    y_pull = y_res / np.sqrt(updated_cov)

    if plot_all:
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updated_cov:\n{updated_cov}")
        print(f"pulls: {y_pull}")
        print("\n")

        max_horizontal = max(detector_layers) + 1
        max_vertical = 40

        fig, ax = plt.subplots()

        # Detector
        for d in range(len(detector_layers)):
            ax.plot(
                [detector_layers[d], detector_layers[d]],
                [-max_vertical, max_vertical],
                "g-",
            )
            ax.plot(detector_layers[d], measurments[d], "gx")

        # Trajectories
        c2u.add_traj_to_plot(
            ax,
            start_params,
            max_horizontal,
            propagators.straight_line_propagator_1D,
            "r",
            "Start Trajectory",
            "-",
        )
        c2u.add_traj_to_plot(
            ax,
            updated_params,
            max_horizontal,
            propagators.straight_line_propagator_1D,
            "b",
            "Final Trajectory",
            "-",
        )
        c2u.add_traj_to_plot(
            ax,
            true_params,
            max_horizontal,
            propagators.straight_line_propagator_1D,
            "k",
            "Unsmeared True Trajectory",
            "-.",
        )

        ax.set(xlabel="x", ylabel="y", title="2D-Fit [y,k]")
        ax.legend()

        # fig.savefig("test.png")
        plt.show()

    ## root fit
    params_res_root, params_pulls_root = c2u.root_fit(
        detector_layers, measurments, cov, "[0]", [true_params]
    )

    return y_pull, y_res, updated_cov, params_pulls_root[0], chi2sum


layers = [int(i ** 1.5) for i in range(2, 10)]
mu_p = []
std_p = []
mu_r = []
std_r = []
mu_c = []
std_c = []
for l in layers:
    print(f"layer {l}")
    draws = 10000
    bins = int(np.sqrt(draws))
    y_pul = []
    y_res = []
    y_cov = []
    y_pul_root = []
    chi2sum = []
    for d in range(draws):
        print("")  # set this line when using spyder, to make root work correctly
        y_p, y_r, y_c, y_p_root, c2s = get_pulls(d < 0, l, 0.1)
        y_pul.append(y_p)
        y_res.append(abs(y_r))
        y_cov.append(np.sqrt(y_c))
        y_pul_root.append(y_p_root)
        chi2sum.append(c2s)

    mu, std = norm.fit(y_pul)
    mu_p.append(mu)
    std_p.append(std)

    mu, std = norm.fit(y_res)
    mu_r.append(mu)
    std_r.append(std)

    mu, std = norm.fit(y_cov)
    mu_c.append(mu)
    std_c.append(std)

    if l == layers[-1] or True:
        c2u.plot_pull_distribution(y_pul, f"y_pulls ({l} hits)")
        c2u.plot_pull_distribution(y_pul_root, f"y_pulls_root ({l} hits)")
        c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y], {l} hits)")

fig, ax = plt.subplots()
ax.plot(layers, std_p, "x")
x_fit = np.linspace(0, max(layers) * 1.2, 100)
ax.set(
    xlabel="n_hits",
    ylabel="std of pulls",
    title="pull distribution depending on $n_{hits}$",
)
ax.legend()

fig.savefig("std_of_pulls.png")
plt.show()

fig, ax = plt.subplots()
ax.plot(layers, mu_r, "x")
params_poly_fit = np.polyfit(1 / np.sqrt(layers), mu_r, 1)
x_fit = np.linspace(2, max(layers) * 1.2, 100)
ax.plot(
    x_fit,
    1 / np.sqrt(x_fit) * params_poly_fit[0] + params_poly_fit[1],
    label="fit: $1/\sqrt{n_{hits}}$",
)
ax.set(
    xlabel="n_hits",
    ylabel="mean of residuals",
    title="residuals depending on $n_{hits}$",
)
ax.legend()

fig.savefig("res.pdf")
plt.show()


fig, ax = plt.subplots()
ax.plot(layers, mu_c, "x")
params_poly_fit = np.polyfit(1 / np.sqrt(layers), mu_c, 1)
x_fit = np.linspace(2, max(layers) * 1.2, 100)
ax.plot(
    x_fit,
    1 / np.sqrt(x_fit) * params_poly_fit[0] + params_poly_fit[1],
    label="fit: $1/\sqrt{n_{hits}}$",
)
ax.set(
    xlabel="n_hits",
    ylabel="mean of $\sigma$",
    title="$\sigma$ depending on $n_{hits}$",
)
ax.legend()

fig.savefig("std.pdf")
plt.show()
