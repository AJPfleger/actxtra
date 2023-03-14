# Alexander J. Pfleger
# 2023-02-28
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
from scipy.stats import norm
import ROOT

import chi2_utilities as c2u


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

    return y_vec#[geo_scatter_sigma == 0]

def get_pulls(plot_all, layers=12, cov=0.1, scatter_sigma_rad=0.05):

    # Parameter
    start_params = np.array([0.0, 3])  # [y, k]
    delta_params = np.zeros_like(start_params)

    # Geometry
    geo_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5, 7, 10, 12, 14, 14.5, 17, 19, 20])
    # geo_layers = np.random.uniform(1, 12, layers)
    # geo_layers.sort()
    geo_scatter_sigma = np.zeros_like(geo_layers)
    geo_scatter_sigma[5] = scatter_sigma_rad
    geo_scatter_sigma[9] = scatter_sigma_rad
    
    detector_layers = geo_layers[geo_scatter_sigma == 0]

    # Hits
    true_params = [12.345, 3]
    measurments_all, cov_meas, measurments_raw = c2u.generate_hits_scatter(
        geo_layers, geo_scatter_sigma, true_params, straight_line_propagator_stepwise, cov, True
    )
    measurments = measurments_all[geo_scatter_sigma == 0]

    updated_params = start_params

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

    return y_pull_root, k_pull_root


scatter_sigma_rad_vec = np.linspace(1e-5,1e-2,100)
draws = 1000
layers = 12
bins = int(np.sqrt(draws))
y_pull_mean = []
y_pull_std = []
k_pull_mean = []
k_pull_std = []

for scatter_sigma_rad in scatter_sigma_rad_vec:
    y_pul_root = []
    k_pul_root = []
    for d in range(draws):
        print("") # set this line when using spyder, to make root work correctly
        y_p_root, k_p_root = get_pulls(False, layers, 0.1, scatter_sigma_rad)
        y_pul_root.append(y_p_root)
        k_pul_root.append(k_p_root)
    
    mu, std = norm.fit(y_pul_root)
    y_pull_mean = np.append(y_pull_mean, mu)
    y_pull_std = np.append(y_pull_std, std)
    mu, std = norm.fit(k_pul_root)
    k_pull_mean = np.append(k_pull_mean, mu)
    k_pull_std = np.append(k_pull_std, std)

fig, ax = plt.subplots()
# plt.hist(chi2sum, bins=bins, density=True)
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 201)
# p = chi2.pdf(x, df, loc, scale)
plt.plot(scatter_sigma_rad_vec, y_pull_std, "b.", label="y-pulls $\sigma$")
plt.plot(scatter_sigma_rad_vec, y_pull_mean, "c.", label="y-pulls $\mu$")
plt.plot(scatter_sigma_rad_vec, k_pull_std, "r.", label="k-pulls $\sigma$")
plt.plot(scatter_sigma_rad_vec, k_pull_mean, "m.", label="k-pulls $\mu$")
# plt.title(f"{title}: k = {df:.3f}, loc = {loc:.3f}, scale = {scale:.3f}")
ax.set(
    # xscale="log",
    xlabel="$\sigma$ of scatter angle / rad",
    ylabel="pulls of each 1e3 draws",
    title="straight line fit with 2 layers of scattering material",
)
ax.legend()
fig.savefig("yk-2_scattering-straight_fit-pull_std_vs_scatter_sigma_lin.pdf")
plt.show()

# c2u.plot_pull_distribution(y_pul_root, f"y_pulls_root ({layers} hits)")
# c2u.plot_pull_distribution(k_pul_root, f"k_pulls_root ({layers} hits)")

