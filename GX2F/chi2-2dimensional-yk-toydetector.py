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
# from scipy.stats import norm
from scipy.optimize import curve_fit

import chi2_utilities as c2u


def fit_func(x, y0, k0):
    return y0 + k0 * x


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
    y_pull = y_res / np.sqrt(updated_cov[0][0])
    k_pull = k_res / np.sqrt(updated_cov[1][1])

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

    popt, pcov = curve_fit(fit_func, detector_layers, measurments)
    y_pull_ref = (popt[0]-true_params[0])/np.sqrt(pcov[0][0])
    k_pull_ref = (popt[1]-true_params[1])/np.sqrt(pcov[1][1])

    return y_pull, k_pull, y_res, k_res, y_pull_ref, k_pull_ref, chi2sum


draws = 1000
layers = 12
bins = int(np.sqrt(draws))
y_pul = []
k_pul = []
y_res_vec = []
k_res_vec = []
y_pul_ref = []
k_pul_ref = []
chi2sum = []
for d in range(draws):
    y_p, k_p, y_res, k_res, y_p_ref, k_p_ref, c2s = get_pulls(d < 5, layers)
    y_pul.append(y_p)
    k_pul.append(k_p)
    y_res_vec.append(y_res)
    k_res_vec.append(k_res)
    y_pul_ref.append(y_p_ref)
    k_pul_ref.append(k_p_ref)
    chi2sum.append(c2s)


c2u.plot_pull_distribution(y_pul, f"y_pulls ({layers} hits)")
c2u.plot_pull_distribution(k_pul, f"k_pulls ({layers} hits)")
c2u.plot_pull_distribution(y_pul_ref, f"y_pulls_reference ({layers} hits)")
c2u.plot_pull_distribution(k_pul_ref, f"k_pulls_reference ({layers} hits)")
c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y,k], {layers} hits)")


from matplotlib.patches import Ellipse

cov = np.cov(y_res_vec, k_res_vec)
lambda_, v = np.linalg.eig(cov)
lambda_ = np.sqrt(lambda_)

fig, ax = plt.subplots()
n_scatter_points = 1000
ax.scatter(y_res_vec[:n_scatter_points], k_res_vec[:n_scatter_points])

for j in range(1, 4):
    ell = Ellipse(xy=(np.mean(y_res_vec), np.mean(k_res_vec)),
                  width=lambda_[0]*j*2, height=lambda_[1]*j*2,
                  angle=np.rad2deg(np.arctan2(*v[:,0][::-1])))
    ell.set_facecolor('none')
    ell.set_edgecolor("red")
    ax.add_artist(ell)

# ax.set_title("y_res vs. k_res")
ax.set(xlabel="y_res", ylabel="k_res", title="y_res vs. k_res, 1e3 runs")
# fig.savefig("yk-y_res-vs-k_res.pdf")
plt.show()
