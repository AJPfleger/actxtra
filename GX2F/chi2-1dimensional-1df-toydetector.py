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


def chi2_1D(V, r):
    return r * (1 / V) * r


def straight_line_propagator(params, x_vec):
    y_vec = np.ones_like(x_vec) * params
    return y_vec


def generate_hits(geometry, true_params, cov=0.1):
    measurments_raw = straight_line_propagator(true_params, geometry)

    cov_meas = [cov] * len(measurments_raw)

    measurments = []
    for mi in range(len(measurments_raw)):
        m = np.random.normal(measurments_raw[mi], np.sqrt(cov_meas[mi]))
        # m = measurments_raw[mi]
        measurments.append(m)

    return measurments, cov_meas, measurments_raw


def get_pulls(plot_all, layers=12, cov=0.1):

    ## Initialising
    n_update = 2

    # Parameter
    start_params = 0.0  # [y0]
    delta_params = 0

    # Geometry
    # detectorLayers = np.array([2, 3, 5, 5.5, 5.7, 6.5,10,10.1,50,98,99,100])
    detectorLayers = np.random.uniform(1, 100, layers)

    # Hits
    true_params = 12.345
    # true_params = [np.random.uniform(-9,9)]
    measurments, cov_meas, measurments_raw = generate_hits(
        detectorLayers, true_params, cov
    )

    updated_params = start_params

    if plot_all:
        maxHorizontal = max(detectorLayers) + 1
        maxVertical = 40
        ## Set up plotting
        fig, ax = plt.subplots()

        def add_traj_to_plot(
            params, color="b", label_text="", style="-", maxHorizontal=maxHorizontal
        ):
            traj = np.array(
                [
                    [0, maxHorizontal],
                    straight_line_propagator(params, [0, maxHorizontal]),
                ]
            )

            ax.plot(0, params, "x" + color)
            ax.plot(traj[0], traj[1], color + style, label=label_text)

    ## Iterating and updating parameters
    for _ in range(n_update):

        updated_params = updated_params + delta_params

        # Iterate over surfaces
        a = 0
        b = 0
        chi2sum = 0
        for d in range(len(detectorLayers)):
            h = detectorLayers[d]
            Vi = cov_meas[d]
            ri = measurments[d] - straight_line_propagator(updated_params, h)
            chi2sum += chi2_1D(Vi, ri)

            ai = 1 / Vi
            bi = ri / Vi

            a += ai
            b += bi

        delta_params = b / a

    updatedCov = 1 / a

    y_res = updated_params - true_params
    y_cov = updatedCov
    y_pull = (updated_params - true_params) / np.sqrt(updatedCov)

    if plot_all:
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updatedCov:\n{updatedCov}")
        print(f"pulls: {y_pull}")
        print("\n")

        # continue plotting
        # Detector
        for d in range(len(detectorLayers)):
            ax.plot(
                [detectorLayers[d], detectorLayers[d]],
                [-maxVertical, maxVertical],
                "g-",
            )
            ax.plot(detectorLayers[d], measurments[d], "gx")

        # Trajectoris
        add_traj_to_plot(start_params, "r", "Start Trajectory", "-")
        add_traj_to_plot(updated_params, "b", "Final Trajectory", "-")
        add_traj_to_plot(true_params, "k", "Unsmeared True Trajectory", "-.")

        ax.set(xlabel="horizontal", ylabel="x", title="1D-Fit")
        ax.legend()

        fig.savefig("setup.pdf")
        plt.show()

    return y_pull, y_res, y_cov


layers = [int(i ** 1.5) for i in range(2, 10)]
mu_p = []
std_p = []
mu_r = []
std_r = []
mu_c = []
std_c = []
for l in layers:
    print(f"layer {l}")
    draws = 1000
    bins = int(np.sqrt(draws))
    y_pul = []
    y_res = []
    y_cov = []
    for d in range(draws):
        y_p, y_r, y_c = get_pulls(d < 0, l, 0.1)
        y_pul.append(y_p)
        y_res.append(abs(y_r))
        y_cov.append(np.sqrt(y_c))

    mu, std = norm.fit(y_pul)
    mu_p.append(mu)
    std_p.append(std)

    mu, std = norm.fit(y_res)
    mu_r.append(mu)
    std_r.append(std)

    mu, std = norm.fit(y_cov)
    mu_c.append(mu)
    std_c.append(std)


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


# plt.hist(y_pulls, bins=bins, density=True)
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 100)
# p = norm.pdf(x, mu, std)
# plt.plot(x, p, "k")
# plt.title(f"y_pulls: mu = {mu:.5f}, std = {std:.5f}")
# plt.show()


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
