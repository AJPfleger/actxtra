# Alexander J. Pfleger
# 2023-03-20
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case
# params are [y, phi]
# scattering is implemented via active surfaces, that give a gaussian smearing
# to the direction. This example includes 2 scattering surfaces

import numpy as np
import matplotlib.pyplot as plt

import chi2_utilities as c2u
import propagators


def get_pulls(plot_all, layers=12, cov=0.1, scatter_sigma_rad=0.0316):

    ## Initialising
    n_update = 250

    # Parameter
    start_params = np.array([0.0, 0.0, 0.0])  # [y, k, theta1]
    # start_params = np.array([12.345, 1, 0.0, 0.0])  # [y, k, theta1, theta2]
    delta_params = np.zeros_like(start_params)

    # Geometry
    geo_layers = np.array(
        [2, 3, 5, 5.5, 5.7, 6.5, 7, 10, 12, 14, 14.01, 17, 19, 20, 21, 22, 23, 24, 25]
    )
    # geo_layers = np.random.uniform(1, 12, layers)
    # geo_layers.sort()
    cov_meas = np.ones_like(geo_layers) * cov
    geo_scatter_sigma = np.zeros_like(geo_layers)
    geo_scatter_sigma[[10]] = scatter_sigma_rad
    # geo_scatter_sigma[[5, 9]] = scatter_sigma_rad

    # detector_layers = geo_layers[geo_scatter_sigma == 0]
    # scatter_layers = geo_layers[geo_scatter_sigma != 0]
    # print(scatter_layers)

    # Hits
    true_params = [0, 0.79]
    scatter_params = np.array([], dtype=float)
    for s_sig in geo_scatter_sigma:
        if s_sig:
            theta = scatter(s_sig)
            # theta = 0.2
            scatter_params = np.append(scatter_params, theta)

    true_params = np.append(true_params, scatter_params)
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments_all, measurments_raw = c2u.generate_hits_scatter(
        geo_layers,
        geo_scatter_sigma,
        true_params,
        propagators.straight_line_propagator_stepwise_2D_scatter_yphi,
        cov_meas,
    )
    # measurments = measurments_all[geo_scatter_sigma == 0]
    # measurments_all = measurments_raw
    updated_params = start_params

    ## Iterating and updating parameters
    for n in range(n_update):

        updated_params = updated_params + delta_params
        
        updated_params[1] = c2u.map_angle_to_right_half(updated_params[1], 0)
        updated_params[2] = c2u.map_angle_to_right_half(
            updated_params[2], updated_params[1]
        )

        predicted_hits = propagators.straight_line_propagator_stepwise_2D_scatter_yphi(
            updated_params, geo_layers, geo_scatter_sigma
        )

        # Iterate over surfaces
        a = np.zeros([len(start_params), len(start_params)])
        b = np.zeros_like(start_params)
        chi2sum = 0
        i_s = 0
        # theta_sum = 0
        # yn_shift = 0
        phi = updated_params[1]
        theta = updated_params[2]
        cosphi2 = np.cos(phi) ** 2
        cosphitheta2 = np.cos(phi + theta) ** 2

        # tant1 = np.tan(t1)
        x_s = np.array([0, 0])
        # dkidt = np.array([0,0])
        for g in range(len(geo_layers)):
            x = geo_layers[g]
            # print(f"x_s = {x_s}")

            if geo_scatter_sigma[g]:  # Scatter layer
                x_s[i_s] = x

                ai = np.zeros([len(start_params), len(start_params)])
                ai[2 + i_s, 2 + i_s] = 1 / geo_scatter_sigma[g] ** 2

                bi = np.zeros_like(start_params)
                bi[2 + i_s] = updated_params[2 + i_s] / geo_scatter_sigma[g] ** 2

                a += ai
                b += bi

                # yn_shift += x * df_dk(k0, theta_sum)
                # theta_sum += updated_params[2+i_s]
                # yn_shift -= x * df_dk(k0, theta_sum)

                # dkidt[i_s] = df_dt(k0,theta_sum)

                chi2sum += updated_params[2 + i_s] ** 2 / geo_scatter_sigma[g] ** 2

                i_s += 1
            else:  # Detector layer
                Vi = cov_meas[g]
                ri = measurments_all[g] - predicted_hits[g]
                chi2sum += c2u.chi2_1D(Vi, ri)

                if i_s == 0:
                    dydy0 = 1
                    dydp0 = x / cosphi2
                    dydt1 = 0
                elif i_s == 1:
                    dydy0 = 1
                    dydp0 = x_s[0] / cosphi2 + (x - x_s[0]) / cosphitheta2
                    dydt1 = (x - x_s[0]) / cosphitheta2
                else:
                    print(f"i_s = {i_s} should not happen")
                # print(dydt1)
                abi_vec = np.array([[dydy0, dydp0, dydt1]])
                # abi_vec = np.array([[dydy0, dydk0, dydt1, dydt2]])

                ai = np.matmul(abi_vec.T, abi_vec)

                # # Second derivatives
                # # It seems that the don't change the result. Only the convergence speed changes.
                # # It seems to be around 20% faster on average, but there is a larger spread in
                # # the convergence time ditribution.
                # if i_s == 0:
                #     # ai[0,0] -= 0
                #     # ai[0,1] -= 0
                #     # ai[1,0] -=
                #     ai[1,1] -= 2 * x * np.tan(phi) / cosphi2 * ri
                #     # ai[0,2] -= 0
                #     # ai[2,0] -= 0
                #     # ai[1,2] -= 0
                #     # ai[2,1] -= 0
                #     # ai[2,2] -= 0
                # elif i_s == 1:
                #     # ai[0,0] -= 0
                #     # ai[0,1] -= 0
                #     # ai[1,0] -= 0
                #     ai[1,1] -= 2 * x_s[0] * np.tan(phi) / cosphi2 * ri + 2 * (x - x_s[0]) * np.tan(phi+theta) / cosphitheta2 * ri
                #     # ai[0,2] -= 0
                #     # ai[2,0] -= 0
                #     ai[1,2] -= 2 * (x - x_s[0]) * np.tan(phi+theta) / cosphitheta2 * ri
                #     ai[2,1] -= 2 * (x - x_s[0]) * np.tan(phi+theta) / cosphitheta2 * ri
                #     ai[2,2] -= 2 * (x - x_s[0]) * np.tan(phi+theta) / cosphitheta2 * ri
                # else:
                #     print(f"i_s = {i_s} should not happen")

                ai /= Vi
                bi = ri / Vi * abi_vec[0]

                a += ai
                b += bi
        # print(f"theta_sum:\n{theta_sum}")
        delta_params = np.linalg.solve(a, b.transpose())

        if abs(delta_params).sum() < 1e-4:
            # print(f"\nmax updates = {n}")
            break

    updated_params[1] = c2u.map_angle_to_right_half(updated_params[1], 0)
    updated_params[2] = c2u.map_angle_to_right_half(
        updated_params[2], updated_params[1]
    )

    updated_cov = np.linalg.inv(a)
    params_res = updated_params - true_params
    y_res, k_res, theta1 = params_res
    # y_res, k_res, theta1, theta2 = params_res
    y_pull = y_res / np.sqrt(updated_cov[0][0])
    k_pull = k_res / np.sqrt(updated_cov[1][1])
    params_pulls = np.zeros_like(params_res)

    for p in range(len(params_res)):
        params_pulls[p] = params_res[p] / np.sqrt(updated_cov[p][p])

    # if params_pulls[2] > 40:
    #     print(f"found\n{a}\n")

    if plot_all or params_res[2] < -1:  # params_pulls[2] > 10 or n > 90 or
        print(f"\nmax updates = {n}")
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
        ax.plot(
            np.append(0, geo_layers),
            np.append(updated_params[0], predicted_hits),
            "b-",
            label="Unsmeared True Trajectory",
        )
        ax.plot(
            np.append(0, geo_layers),
            np.append(true_params[0], measurments_raw),
            "k-",
            label="Unsmeared True Trajectory",
        )

        ax.set(xlabel="x", ylabel="y", title="2D-Fit [y,k]")
        ax.legend()

        # fig.savefig("toydetector-scattering-straight-fit.pdf")
        plt.show()
        print(f"delta_params = {delta_params}")
    return params_res, params_pulls, n  # chi2sum


np.random.seed(10)
draws = 10000
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
    # print(f"draw {d}")
    # print("") # set this line when using spyder, to make root work correctly
    p_res, p_pulls, c2s = get_pulls(d < 5, layers)
    # if d == 645:
    #     continue
    # if not -0.1 < p_res[2] < 0.1:
    #     continue
    y_pul.append(p_pulls[0])
    k_pul.append(p_pulls[1])
    # if p_pulls[2] < 50:
    t1_pul.append(p_pulls[2])
    # t2_pul.append(p_pulls[3])
    y_res.append(p_res[0])
    k_res.append(p_res[1])
    t1_res.append(p_res[2])
    # t2_res.append(p_res[3])
    chi2sum.append(c2s)


c2u.plot_pull_distribution(y_res, f"y_res ({layers} hits)")
c2u.plot_pull_distribution(k_res, f"phi_res({layers} hits)")
c2u.plot_pull_distribution(t1_res, f"t1_res ({layers} hits)")
# c2u.plot_pull_distribution(t2_res, f"t2_res ({layers} hits)")

c2u.plot_pull_distribution(y_pul, f"y_pulls ({layers} hits)")
c2u.plot_pull_distribution(k_pul, f"phi_pulls ({layers} hits)")
c2u.plot_pull_distribution(t1_pul, f"t1_pulls ({layers} hits)")
# c2u.plot_pull_distribution(t2_pul, f"t2_pulls ({layers} hits)")
c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y,k], {layers} hits)")

# plt.plot(t1_res,".")
# plt.ylim(-0.5, 0.2)


c2u.plot_pull_distribution(chi2sum, f"y-phi-1S without second derivatives (18 hits)")

plt.plot(chi2sum, ".")
plt.ylim(0, 180)

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
