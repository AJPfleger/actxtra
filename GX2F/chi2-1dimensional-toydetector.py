# Alexander J. Pfleger
# 2022-01-24
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 1D case with just a single degree of freedom
# params are [y]

import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy.stats import norm

import chi2_utilities as c2u
import propagators


def ai_bi(ri, Vi, xi, params):

    ai = 1 / Vi * np.array([[1]])
    bi = ri / Vi * np.array([1])

    return ai, bi


def get_pulls(plot_all, layers=12, cov=0.1):

    ## Initialising
    n_update = 2
    true_params = [1.2345]
    # true_params = [np.random.uniform(-9,9)]
    start_params = [1.0]  # [y0]

    # Geometry
    # detector_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5,10,10.1,50,98,99,100])
    detector_layers = np.random.uniform(1, 100, layers)

    # Generate hits
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detector_layers, true_params, propagators.straight_line_propagator_1D, cov
    )

    # Fit
    a, updated_params, chi2sum, updated_cov = c2u.gx2f(
        start_params,
        detector_layers,
        cov_meas,
        measurments,
        propagators.straight_line_propagator_1D,
        n_update,
        ai_bi,
    )

    params_res, params_pulls = c2u.calc_res_pulls(
        updated_params, true_params, updated_cov
    )

    y_res = params_res[0]
    y_pull = params_pulls[0]

    if plot_all or True:
        geo_scatter_sigma = detector_layers * 0
        start_traj = propagators.straight_line_propagator_stepwise_2D_scatter_yk(
            [start_params[0], 0], detector_layers, geo_scatter_sigma,
        )
        predicted_hits = propagators.straight_line_propagator_stepwise_2D_scatter_yk(
            [updated_params[0], 0], detector_layers, geo_scatter_sigma,
        )

        c2u.plot_current_state(
            updated_params,
            true_params,
            a,
            updated_cov,
            measurments,  # measurments_all,
            detector_layers,  # geo_layers,
            geo_scatter_sigma,
            predicted_hits,
            measurments_raw,
            "",
            params_pulls,
            "1D-Fit [y]",
            "",
            start_params,
            start_traj,
        )

    ## root fit
    params_res_root, params_pulls_root = c2u.root_fit(
        detector_layers, measurments, cov, "[0]", true_params
    )

    return y_pull, y_res, updated_cov, params_pulls_root[0], chi2sum[0]


logging.getLogger().setLevel(logging.INFO)
np.random.seed(10)
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
