# Alexander J. Pfleger
# 2022-11-25
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case

import numpy as np
import matplotlib.pyplot as plt

# from scipy.optimize import curve_fit
import ROOT

import chi2_utilities as c2u
import propagators


def get_pulls(plot_all, layers=5, cov=0.1, phi_true=-1):

    ## Initialising
    n_update = 15

    # Parameter
    start_params = np.array([0.0, 0])  # [y, phi]
    delta_params = np.zeros_like(start_params)

    # Geometry
    # detector_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5, 10, 50, 98, 99, 100])
    detector_layers = np.random.uniform(2, 10, layers)

    # Hits
    true_params = [12.345, phi_true]
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detector_layers,
        true_params,
        propagators.straight_line_propagator_2D_yphi,
        0.1,
        True,
    )
    # proj = np.array([[1,0]]) # projects onto x
    start_params = true_params
    updated_params = start_params

    ## Iterating and updating parameters
    for _ in range(n_update):

        updated_params = updated_params + delta_params
        updated_params[1] = c2u.map_angle_to_right_half(updated_params[1], 0)

        # Iterate over surfaces
        a = np.zeros([2, 2])
        b = np.zeros_like(start_params)
        chi2sum = 0
        for d in range(len(detector_layers)):
            x = detector_layers[d]
            Vi = cov_meas[d]
            ri = measurments[d] - propagators.straight_line_propagator_2D_yphi(
                updated_params, x
            )
            chi2sum += c2u.chi2_1D(Vi, ri)

            x_cos2 = x / np.cos(updated_params[1]) ** 2
            d2chi_dphi2 = (
                x_cos2 ** 2
            )  # - 2 * ri / h_cos2 ** 2 * np.tan(updated_params[1]))
            ai = 1 / Vi * np.array([[1, x_cos2], [x_cos2, d2chi_dphi2],])
            bi = ri / Vi * np.array([1, x_cos2])

            a += ai
            b += bi

        delta_params = np.linalg.solve(a, b.transpose())

    updated_cov = np.linalg.inv(a)

    params_res, params_pulls = c2u.calc_res_pulls(
        updated_params, true_params, updated_cov
    )

    y_res, phi_res = params_res
    y_pull, phi_pull = params_pulls

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
            propagators.straight_line_propagator_2D_yphi,
            "r",
            "Start Trajectory",
            "-",
        )
        c2u.add_traj_to_plot(
            ax,
            updated_params,
            max_horizontal,
            propagators.straight_line_propagator_2D_yphi,
            "b",
            "Final Trajectory",
            "-",
        )
        c2u.add_traj_to_plot(
            ax,
            true_params,
            max_horizontal,
            propagators.straight_line_propagator_2D_yphi,
            "k",
            "Unsmeared True Trajectory",
            "-.",
        )

        ax.set(xlabel="x", ylabel="y", title="2D-Fit [y,k]")
        ax.legend()

        # fig.savefig("test.png")
        plt.show()

    ## root fit
    x_root = np.array(detector_layers)
    y_root = np.array(measurments)
    ex_root = x_root * 0
    ey_root = ex_root + np.sqrt(cov)
    g_root = ROOT.TGraphErrors(len(x_root), x_root, y_root, ex_root, ey_root)
    f_root = ROOT.TF1("f_root", "[0] + tan([1])*x")
    t_root = g_root.Fit(f_root, "S")
    # print(t_root)

    y_res_root = t_root.Parameter(0) - true_params[0]
    y_std_root = t_root.Error(0)
    y_pull_root = y_res_root / y_std_root

    phi_root = c2u.map_angle_to_right_half(t_root.Parameter(1), 0)

    phi_res_root = phi_root - true_params[1]
    phi_std_root = t_root.Error(1)
    phi_pull_root = phi_res_root / phi_std_root

    chi2sum = abs(1 - y_pull / y_pull_root)
    updated_cov[1][1] = abs(1 - phi_pull / phi_pull_root)
    if abs(y_pull_root - y_pull) > 1e-3:
        print("----------e-r-r-o-r----------")
        print(y_pull_root - y_pull)

    return (
        y_pull,
        phi_pull,
        updated_cov[1][1],
        y_res_root,
        phi_res_root,
        y_pull_root,
        phi_pull_root,
        chi2sum,
    )


draws = 10000
layers = 30
bins = int(np.sqrt(draws))
y_pul = []
phi_pul = []
phi_covs = []
y_res_vec = []
phi_res_vec = []
y_pul_ref = []
phi_pul_ref = []
chi2sum = []
for d in range(draws):
    print("")  # set this line when using spyder, to make root work correctly
    y_p, phi_p, phi_c, y_res, phi_res, y_p_ref, phi_p_ref, c2s = get_pulls(
        d < 5, layers, 0.1, 1
    )
    y_pul.append(y_p)
    phi_pul.append(phi_p)
    phi_covs.append(phi_c)
    y_res_vec.append(y_res)
    phi_res_vec.append(phi_res)
    y_pul_ref.append(y_p_ref)
    phi_pul_ref.append(phi_p_ref)
    chi2sum.append(c2s)

c2u.plot_pull_distribution(y_pul, f"y_pulls ({layers} hits)")
c2u.plot_pull_distribution(y_pul_ref, f"y_pulls_root ({layers} hits)")
c2u.plot_pull_distribution(phi_pul, f"$\phi$_pulls ({layers} hits)")
c2u.plot_pull_distribution(phi_pul_ref, f"$\phi$_pulls_root ({layers} hits)")
c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y], {layers} hits)")
c2u.plot_chi2_distribution(phi_covs, f"$\sigma_\phi$ ([y], {layers} hits)")

from scipy.stats import skewnorm

a, loc, scale = skewnorm.fit(phi_pul)
plt.hist(phi_pul, bins=bins, density=True)
xminmax = 5
x = np.linspace(-xminmax, xminmax, 100)
plt.xlim([-xminmax, xminmax])
p = skewnorm.pdf(x, a, loc, scale)
plt.plot(x, p, "k")
delta = a / np.sqrt(1 + a ** 2)
mu = loc + scale * delta * np.sqrt(2 / np.pi)
std = scale * np.sqrt(1 - 2 * delta ** 2 / np.pi)
skew = (
    (4 - np.pi)
    / 2
    * (delta * np.sqrt(2 / np.pi)) ** 3
    / (1 - 2 * delta ** 2 / np.pi) ** (3 / 2)
)
plt.title(f"phi_pulls: mu = {mu:.5f}, std = {std:.5f}, skew = {skew:.5f}")
plt.show()


from matplotlib.patches import Ellipse

cov = np.cov(y_res_vec, phi_res_vec)
lambda_, v = np.linalg.eig(cov)
lambda_ = np.sqrt(lambda_)

fig, ax = plt.subplots()
n_scatter_points = 10000
ax.scatter(y_res_vec[:n_scatter_points], phi_res_vec[:n_scatter_points])

for j in range(1, 4):
    ell = Ellipse(
        xy=(np.mean(y_res_vec), np.mean(phi_res_vec)),
        width=lambda_[0] * j * 2,
        height=lambda_[1] * j * 2,
        angle=np.rad2deg(np.arctan2(*v[:, 0][::-1])),
    )
    ell.set_facecolor("none")
    ell.set_edgecolor("red")
    ax.add_artist(ell)

ax.set(xlabel="y_res", ylabel="phi_res", title="y_res vs. phi_res, 1e4 runs")
# fig.savefig("yk-y_res-vs-phi_res.pdf")
plt.show()


# # #### vvvv experimental vvvv ####
# from scipy.stats import norm
# from scipy.stats import skewnorm

# skew_vec = []

# draws = 10000
# phi_test_vec = np.linspace(-1,1,100)
# for phi_test in phi_test_vec:
#     bins = int(np.sqrt(draws))
#     # y_pulls = []
#     phi_pulls = []
#     # phi_covs = []
#     # y_res_vec = []
#     # phi_res_vec = []
#     for d in range(draws):
#         print("") # set this line when using spyder, to make root work correctly
#         y_p, phi_p, phi_c, y_res, phi_res, y_p_ref, phi_p_ref, c2s = get_pulls(False, layers=5, cov=0.1, phi_true=phi_test)
#         # y_pulls.append(y_p)
#         phi_pulls.append(phi_p)
#         # phi_covs.append(phi_c)
#         # y_res_vec.append(y_res)
#         # phi_res_vec.append(phi_res)


#     # fig, ax = plt.subplots()
#     a, loc, scale = skewnorm.fit(phi_pulls)
#     # plt.hist(phi_pulls, bins=bins, density=True)
#     xminmax = 5
#     # x = np.linspace(-xminmax, xminmax, 100)
#     plt.xlim([-xminmax,xminmax])
#     # p = skewnorm.pdf(x, a, loc, scale)

#     delta = a/np.sqrt(1 + a**2)
#     mu = loc + scale * delta * np.sqrt(2/np.pi)
#     std = scale * np.sqrt(1 - 2*delta**2/np.pi)
#     skew = (4-np.pi)/2*(delta*np.sqrt(2/np.pi))**3/(1 - 2*delta**2/np.pi)**(3/2)

#     skew_vec.append(skew)

#     # if abs(skew) < 0.01:
#     #     plt.plot(x, p, "r")
#     # else:
#     #     plt.plot(x, p, "k")
#     # ax.set(title=f"phi_true={phi_test}, mu={mu:.4f}, std={std:.4f}, skew={skew:.4f}")


#     # plt.show()
# from scipy.special import erf

# fig, ax = plt.subplots()
# plt.plot(phi_test_vec, skew_vec,".")
# # skew_fit = [0]
# # for i in range(len(skew_vec)-1):
#     # skew_fit.append(skew_fit[i] + skew_vec[i])
# skew_fit = (4-np.pi)/2*erf(phi_test_vec/np.sqrt(2*0.5))
# plt.plot(phi_test_vec, skew_fit,"g-")
# plt.plot([-1,1],[0,0],"r-")
# plt.plot([0,0],[-0.3,0.3],"r-")
# ax.set(xlabel="phi_true", ylabel="skew_phi_pulls", title="phi_true vs. skew_phi_pulls, 1e4 runs, ROOT.TF1")
# # ax.set(xlabel="phi_true", ylabel="skew_phi_pulls", title="phi_true vs. skew_phi_pulls, 1e4 runs, phi_start=true, nUpdates=15")
# fig.savefig("yphi_skew_of_pulls_phitrue_root_fit.pdf")
# plt.show()
