# Alexander J. Pfleger
# 2023-02-23
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case
# params are [y, k]
# scattering is implemented via active surfaces, that give a gaussian smearing
# to the direction. This example includes 2 scattering surfaces

import numpy as np
import matplotlib.pyplot as plt

# import ROOT

import chi2_utilities as c2u


def fit_func(x, y0, k0):
    return y0 + k0 * x


def straight_line_propagator(params, x_vec):
    y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

    return y_vec


# def straight_line_propagator_stepwise(start_params, geo_pos, geo_scatter_sigma):

#     assert(len(geo_pos) == len(geo_scatter_sigma))

#     current_x_y_k = [0, start_params[0], start_params[1]]
#     new_x_y_k = [0, start_params[0], start_params[1]]

#     y_vec = np.zeros_like(geo_pos)
#     theta_vec = np.zeros_like(geo_pos) # Saves the theta at scattering surfaces

#     for g in range(len(geo_pos)):
#         # make hit
#         new_x_y_k = current_x_y_k.copy()
#         new_x_y_k[0] = geo_pos[g]
#         dx = new_x_y_k[0] - current_x_y_k[0]
#         new_x_y_k[1] = current_x_y_k[1] + dx*current_x_y_k[2]
#         y_vec[g] = new_x_y_k[1]

#         if geo_scatter_sigma[g] != 0: # scatter
#             theta_vec[g] = np.arctan(current_x_y_k[2])
#             theta_new = np.random.normal(theta_vec[g], geo_scatter_sigma[g])
#             new_x_y_k[2] = np.tan(theta_new)

#         current_x_y_k = new_x_y_k.copy()

#     # y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

#     return y_vec, theta_vec


def straight_line_propagator_stepwise(params, geo_pos, is_scatter):
    assert len(geo_pos) == len(is_scatter)

    current_x_y_k = [0, params[0], params[1]]
    new_x_y_k = [0, params[0], params[1]]

    scatter_params = params[2:]

    y_vec = np.zeros_like(geo_pos)
    i_s = 0  # to count through the scattering surfaces
    for g in range(len(geo_pos)):
        # make hit
        new_x_y_k = current_x_y_k.copy()
        new_x_y_k[0] = geo_pos[g]

        dx = new_x_y_k[0] - current_x_y_k[0]
        new_x_y_k[1] = current_x_y_k[1] + current_x_y_k[2] * dx

        if is_scatter[g]:
            kn = current_x_y_k[2]
            ks = np.tan(scatter_params[i_s])
            new_x_y_k[2] = (kn + ks) / (1 - kn * ks)
            i_s += 1
        else:
            new_x_y_k[2] = current_x_y_k[2]

        y_vec[g] = new_x_y_k[1]
        current_x_y_k = new_x_y_k.copy()

    return y_vec


def scatter(sigma):
    return np.random.normal(0, sigma)


def df_dk(k0, theta_sum):
    tan = np.tan(theta_sum)
    
    numerator = 1 + tan ** 2
    denominator = (1 - k0 * tan) **2
    
    # (1 + np.tan(theta_sum) ** 2) / (1 - k0 * np.tan(theta_sum)) + yn_shift
    
    return numerator / denominator


def df_dt(k0, theta_sum):
    tan = np.tan(theta_sum)
    cos = np.cos(theta_sum)
    
    numerator = 1 + k0 ** 2
    denominator = (cos * (1 - k0 * tan)) ** 2
    
    return numerator / denominator


def get_pulls(plot_all, layers=12, cov=0.1, scatter_sigma_rad=0.05):

    ## Initialising
    n_update = 15

    # Parameter
    start_params = np.array([0.0, 3, 0.11, 0.11])  # [y, k, theta1, theta2]
    delta_params = np.zeros_like(start_params)

    # Geometry
    geo_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5, 7, 10, 12, 14, 14.5, 17, 19, 20])
    # geo_layers = np.random.uniform(1, 12, layers)
    # geo_layers.sort()
    cov_meas = np.ones_like(geo_layers) * cov
    geo_scatter_sigma = np.zeros_like(geo_layers)
    geo_scatter_sigma[[5, 9]] = scatter_sigma_rad
    
    # detector_layers = geo_layers[geo_scatter_sigma == 0]
    # scatter_layers = geo_layers[geo_scatter_sigma != 0]
    # print(scatter_layers)

    # Hits
    true_params = [12.345, 3]
    scatter_params = np.array([], dtype=float)
    for s_sig in geo_scatter_sigma:
        if s_sig:
            # theta = scatter(s_sig)
            theta = 0.1
            scatter_params = np.append(scatter_params, theta)

    true_params = np.append(true_params, scatter_params)
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments_all, measurments_raw = c2u.generate_hits_scatter(
        geo_layers,
        geo_scatter_sigma,
        true_params,
        straight_line_propagator_stepwise,
        cov_meas,
    )
    # measurments = measurments_all[geo_scatter_sigma == 0]
    # measurments_all = measurments_raw
    updated_params = start_params

    ## Iterating and updating parameters
    for _ in range(n_update):

        updated_params = updated_params + delta_params
        predicted_hits = straight_line_propagator_stepwise(
            updated_params, geo_layers, geo_scatter_sigma
        )

        # Iterate over surfaces
        a = np.zeros([len(start_params), len(start_params)])
        b = np.zeros_like(start_params)
        chi2sum = 0
        i_s = 0
        theta_sum = 0
        yn_shift = 0
        k0 = start_params[1]
        x_s = np.array([0,0])
        dkidt = np.array([0,0])
        for g in range(len(geo_layers)):
            x = geo_layers[g]
            
            if geo_scatter_sigma[g]: # Scatter layer
                x_s[i_s] = x
                
                ai = np.zeros([len(start_params), len(start_params)])
                ai[2+i_s, 2+i_s] = 1 / geo_scatter_sigma[g]**2
                
                bi = np.zeros_like(start_params)
                bi[2+i_s] = updated_params[2+i_s] / geo_scatter_sigma[g]**2
                
                a += ai
                b += bi
                
                yn_shift += x * df_dk(k0, theta_sum)
                theta_sum += updated_params[2+i_s]
                yn_shift -= x * df_dk(k0, theta_sum)
                
                dkidt[i_s] = df_dt(k0,theta_sum)
                
                i_s += 1
            else: # Detector layer
                Vi = cov_meas[g]
                ri = measurments_all[g] - predicted_hits[g]
                chi2sum += c2u.chi2_1D(Vi, ri)
    # missing: y1 terms
                dydy0 = 1
                dydk0 = x * df_dk(k0, theta_sum) + yn_shift

                if i_s == 0:
                    dydt1 = 0
                    dydt2 = 0
                elif i_s == 1:
                    dydt1 = (x - x_s[0]) * dkidt[0]
                    dydt2 = 0
                elif i_s == 2:
                    dydt1 = (x - x_s[1]) * dkidt[1] + (x_s[1] - x_s[0]) * dkidt[0] 
                    dydt2 = (x - x_s[1]) * dkidt[1]
                else:
                    print(f"i_s = {i_s} should not happen")
                print(dydt1)
                abi_vec = np.array([[dydy0, dydk0, dydt1, dydt2]])
                ai = 1 / Vi * np.matmul(abi_vec.T,abi_vec)
                bi = ri / Vi * abi_vec[0]
    
                a += ai
                b += bi
        # print(f"theta_sum:\n{theta_sum}")
        delta_params = np.linalg.solve(a, b.transpose())

    updated_cov = np.linalg.inv(a)
    params_res = updated_params - true_params
    y_res, k_res, theta1, theta2 = params_res
    y_pull = y_res / np.sqrt(updated_cov[0][0])
    k_pull = k_res / np.sqrt(updated_cov[1][1])
    params_pulls = np.zeros_like(params_res)
    for p in range(len(params_res)):
        params_pulls[p] = params_res[p] / np.sqrt(updated_cov[p][p])

    if plot_all:
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updated_cov:\n{updated_cov}")
        print(f"pulls: {y_pull}, {k_pull}")
        print(params_pulls)
        print("\n")

        # max_horizontal = max(geo_layers) + 1
        delta_measurments = abs(max(measurments_all) - min(measurments_all))
        min_vertical = min(measurments_all) - 0.3 * delta_measurments
        max_vertical = max(measurments_all) + 0.3 * delta_measurments

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
        
        ax.plot(geo_layers, measurments_all, "gx")

        # Trajectoris
        # c2u.add_traj_to_plot(ax, start_params, max_horizontal, straight_line_propagator, "r", "Start Trajectory", "-")
        # c2u.add_traj_to_plot(ax, updated_params, max_horizontal, straight_line_propagator, "b", "Final Trajectory", "-")
        ax.plot(np.append(0, geo_layers), np.append(updated_params[0], predicted_hits), "b-", label="Unsmeared True Trajectory")
        ax.plot(np.append(0, geo_layers), np.append(true_params[0], measurments_raw), "k-", label="Unsmeared True Trajectory")

        ax.set(xlabel="x", ylabel="y", title="2D-Fit [y,k]")
        ax.legend()

        #fig.savefig("toydetector-scattering-straight-fit.pdf")
        plt.show()

    return params_res, params_pulls, chi2sum


draws = 1
layers = 12
bins = int(np.sqrt(draws))
y_pul = []
k_pul = []
t1_pul = []
t2_pul = []
y_res = []
k_res = []
t1_res = []
t2_res = []
chi2sum = []
for d in range(draws):
    # print("") # set this line when using spyder, to make root work correctly
    p_res, p_pulls, c2s = get_pulls(d < 5, layers)
    y_pul.append(p_pulls[0])
    k_pul.append(p_pulls[1])
    t1_pul.append(p_pulls[2])
    t2_pul.append(p_pulls[3])
    y_res.append(p_res[0])
    k_res.append(p_res[1])
    t1_res.append(p_res[2])
    t2_res.append(p_res[3])
    chi2sum.append(c2s)

# c2u.plot_pull_distribution(y_res, f"y_res ({layers} hits)")
# c2u.plot_pull_distribution(k_res, f"k_res({layers} hits)")
# c2u.plot_pull_distribution(t1_res, f"t1_res ({layers} hits)")
# c2u.plot_pull_distribution(t2_res, f"t2_res ({layers} hits)")

# c2u.plot_pull_distribution(y_pul, f"y_pulls ({layers} hits)")
# c2u.plot_pull_distribution(k_pul, f"k_pulls ({layers} hits)")
# c2u.plot_pull_distribution(t1_pul, f"t1_pulls ({layers} hits)")
# c2u.plot_pull_distribution(t2_pul, f"t2_pulls ({layers} hits)")
# c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y,k], {layers} hits)")


# from matplotlib.patches import Ellipse

# cov = np.cov(y_res_vec, k_res_vec)
# lambda_, v = np.linalg.eig(cov)
# lambda_ = np.sqrt(lambda_)

# fig, ax = plt.subplots()
# n_scatter_points = 1000
# ax.scatter(y_res_vec[:n_scatter_points], k_res_vec[:n_scatter_points])

# for j in range(1, 4):
#     ell = Ellipse(xy=(np.mean(y_res_vec), np.mean(k_res_vec)),
#                   width=lambda_[0]*j*2, height=lambda_[1]*j*2,
#                   angle=np.rad2deg(np.arctan2(*v[:,0][::-1])))
#     ell.set_facecolor('none')
#     ell.set_edgecolor("red")
#     ax.add_artist(ell)

# ax.set(xlabel="y_res", ylabel="k_res", title="y_res vs. k_res, 1e3 runs")
# # fig.savefig("yk-y_res-vs-k_res.pdf")
# plt.show()
