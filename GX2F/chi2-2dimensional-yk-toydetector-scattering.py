# Alexander J. Pfleger
# 2023-02-23
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case
# params are [y, k]
# scattering is implemented via active surfaces, that give a gaussian smearing
# to the direction

import numpy as np
import matplotlib.pyplot as plt
import ROOT

import chi2_utilities as c2u


def fit_func(x, y0, k0):
    return y0 + k0 * x


def straight_line_propagator(params, x_vec):
    y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

    return y_vec

def straight_line_propagator_stepwise(start_params, geo_pos, geo_scatter_sigma):
    
    assert(len(geo_pos) == len(geo_scatter_sigma))
    
    current_x_y_k = [0, start_params[0], start_params[1]]
    new_x_y_k = [0, start_params[0], start_params[1]]
    
    y_vec = np.zeros_like(geo_pos)
    
    for g in range(len(geo_pos)):
        # make hit
        new_x_y_k = current_x_y_k.copy()
        new_x_y_k[0] = geo_pos[g]
        dx = new_x_y_k[0] - current_x_y_k[0]
        new_x_y_k[1] = current_x_y_k[1] + dx*current_x_y_k[2]
        y_vec[g] = new_x_y_k[1]
        
        if geo_scatter_sigma[g] != 0: # scatter
            phi = np.arctan(current_x_y_k[2])
            phi_new = np.random.normal(phi, geo_scatter_sigma[g])
            new_x_y_k[2] = np.tan(phi_new)
        
        current_x_y_k = new_x_y_k.copy()
    
    # y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

    return y_vec#[geo_scatter_sigma == 0]

def get_pulls(plot_all, layers=12, cov=0.1, scatter_sigma_rad=0.05):

    ## Initialising
    n_update = 15

    # Parameter
    start_params = np.array([0.0, 3])  # [y, k]
    delta_params = np.zeros_like(start_params)

    # Geometry
    geo_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5, 7, 10, 12, 14, 14.5, 17, 19, 20])
    # geo_layers = np.random.uniform(1, 12, layers)
    geo_layers.sort()
    geo_scatter_sigma = np.zeros_like(geo_layers)
    geo_scatter_sigma[5] = scatter_sigma_rad
    geo_scatter_sigma[9] = scatter_sigma_rad
    
    detector_layers = geo_layers[geo_scatter_sigma == 0]

    # Hits
    true_params = [12.345, 3]
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments_all, cov_meas, measurments_raw = c2u.generate_hits_scatter(
        geo_layers, geo_scatter_sigma, true_params, straight_line_propagator_stepwise, cov, True
    )
    measurments = measurments_all[geo_scatter_sigma == 0]

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
        delta_measurments = abs(max(measurments) - min(measurments))
        min_vertical = min(measurments) - 0.3 * delta_measurments
        max_vertical = max(measurments) + 0.3 * delta_measurments

        fig, ax = plt.subplots()

        # Detector
        for d in range(len(geo_layers)):
            if geo_scatter_sigma[d] == 0:
                line_style_surface = "g-"
            else:
                line_style_surface = "r:"
            ax.plot(
                [geo_layers[d], geo_layers[d]],
                [min_vertical, max_vertical],
                line_style_surface,
            )
        
        ax.plot(detector_layers, measurments, "gx")

        # Trajectoris
        # c2u.add_traj_to_plot(ax, start_params, max_horizontal, straight_line_propagator, "r", "Start Trajectory", "-")
        c2u.add_traj_to_plot(ax, updated_params, max_horizontal, straight_line_propagator, "b", "Final Trajectory", "-")
        ax.plot(np.append(0, geo_layers), np.append(true_params[0], measurments_raw), "k-", label="Unsmeared True Trajectory")

        ax.set(xlabel="horizontal", ylabel="x", title="2D-Fit [y,k]")
        ax.legend()

        fig.savefig("toydetector-scattering-straight-fit.pdf")
        plt.show()

    ## root fit
    x_root = np.array(detector_layers)
    y_root = np.array(measurments)
    ex_root = x_root*0
    ey_root = ex_root + np.sqrt(cov)
    g_root = ROOT.TGraphErrors(len(x_root), x_root, y_root, ex_root, ey_root)
    f_root = ROOT.TF1("f_root", "[0] + [1]*x")
    t_root = g_root.Fit(f_root,"S")

    y_res_root = t_root.Parameter(0)-true_params[0]
    y_std_root = t_root.Error(0)
    y_pull_root = y_res_root / y_std_root

    k_res_root = t_root.Parameter(1)-true_params[1]
    k_std_root = t_root.Error(1)
    k_pull_root = k_res_root / k_std_root

    return y_pull, k_pull, y_res, k_res, y_pull_root, k_pull_root, chi2sum


draws = 100
layers = 12
bins = int(np.sqrt(draws))
y_pul = []
k_pul = []
y_res_vec = []
k_res_vec = []
y_pul_root = []
k_pul_root = []
chi2sum = []
for d in range(draws):
    print("") # set this line when using spyder, to make root work correctly
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
    ell = Ellipse(xy=(np.mean(y_res_vec), np.mean(k_res_vec)),
                  width=lambda_[0]*j*2, height=lambda_[1]*j*2,
                  angle=np.rad2deg(np.arctan2(*v[:,0][::-1])))
    ell.set_facecolor('none')
    ell.set_edgecolor("red")
    ax.add_artist(ell)

ax.set(xlabel="y_res", ylabel="k_res", title="y_res vs. k_res, 1e3 runs")
# fig.savefig("yk-y_res-vs-k_res.pdf")
plt.show()
