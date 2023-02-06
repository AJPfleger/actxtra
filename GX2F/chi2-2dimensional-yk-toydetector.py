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
from scipy.stats import norm

import chi2_utilities as c2u


def straight_line_propagator(params, x_vec):
    y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

    return y_vec


def get_pulls(plot_all, layers=12, cov=0.1):

    ## Initialising
    n_update = 15

    # Parameter
    start_params = np.array([0.0, 3])  # [y, k]
    delta_params = np.zeros_like(start_params)

    # Geometry
    # detector_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5])  # ,10,50,98,99,100])
    detector_layers = np.random.uniform(1, 12, layers)

    # Hits
    true_params = [12.345, 1]
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detector_layers, true_params, straight_line_propagator, cov, True
    )

    updated_params = start_params

    ## Iterating and updating parameters
    for _ in range(n_update):

        updated_params = updated_params + delta_params

        # Iterate over surfaces
        a = np.zeros([2, 2])
        b = np.zeros_like(start_params)
        chi2sum = 0
        for d in range(len(detector_layers)):
            h = detector_layers[d]
            # propagatedCov = straightLineCovarianceTransport(updated_cov,updated_params,h)
            Vi = cov_meas[d]
            ri = measurments[d] - straight_line_propagator(updated_params, h)
            chi2sum += c2u.chi2_1D(Vi, ri)

            ai = 1 / Vi * np.array([[1, h], [h, h ** 2],])
            bi = ri / Vi * np.array([1, h])

            a += ai
            b += bi

        delta_params = np.linalg.solve(a, b.transpose())

    updated_cov = np.linalg.inv(a)

    y_res, k_res = updated_params - true_params
    k_res = np.arctan(updated_params[1]) - np.arctan(true_params[1])
    y_pull = y_res / np.sqrt(updated_cov[0][0])
    k_pull = k_res #/ np.sqrt( (updated_cov[1][1]) )

    if plot_all:        
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updated_cov:\n{updated_cov}")
        print(f"pulls: {y_pull}, {k_pull}")
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

        # Trajectoris
        c2u.add_traj_to_plot(ax, start_params, max_horizontal, straight_line_propagator, "r", "Start Trajectory", "-")
        c2u.add_traj_to_plot(ax, updated_params, max_horizontal, straight_line_propagator, "b", "Final Trajectory", "-")
        c2u.add_traj_to_plot(ax, true_params, max_horizontal, straight_line_propagator, "k", "Unsmeared True Trajectory", "-.")

        ax.set(xlabel="horizontal", ylabel="x", title="2D-Fit [y,k]")
        ax.legend()

        # fig.savefig("test.png")
        plt.show()

    return y_pull, k_pull, y_res, k_res


draws = 100000
bins = int(np.sqrt(draws))
y_pulls = []
k_pulls = []
y_res_vec = []
k_res_vec = []
for d in range(draws):
    y_p, k_p, y_res, k_res = get_pulls(d < 5)
    y_pulls.append(y_p)
    k_pulls.append(k_p)
    y_res_vec.append(y_res)
    k_res_vec.append(k_res)


mu, std = norm.fit(y_pulls)
fig, ax = plt.subplots()
ax.hist(y_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
ax.plot(x, p, "k")
ax.set(title=f"y_pulls: mu = {mu:.5f}, std = {std:.5f}")
fig.savefig("yk-y_pulls.png")
plt.show()


mu, std = norm.fit(k_pulls)
fig, ax = plt.subplots()
ax.hist(k_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
ax.plot(x, p, "k")
ax.set(title=f"k_pulls: mu = {mu:.5f}, std = {std:.5f}")
fig.savefig("yk-k_pulls.png")
plt.show()
