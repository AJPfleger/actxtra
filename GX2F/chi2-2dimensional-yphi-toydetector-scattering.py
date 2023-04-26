# Alexander J. Pfleger
# 2023-03-20
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case
# params are [y, phi]
# scattering is implemented via active surfaces, that give a gaussian smearing
# to the direction. This example includes 1 scattering surfaces

import numpy as np
import matplotlib.pyplot as plt
import logging

import chi2_utilities as c2u
import propagators

# scatter-angle for 1 GeV pion with 300 um silicon is about 0.8 mrad
# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwizosSqoLv-AhXUGogKHQYBBKkQFnoECCEQAQ&url=https%3A%2F%2Fe-publishing.cern.ch%2Findex.php%2FCYRSP%2Farticle%2Fdownload%2F534%2F396%2F1738&usg=AOvVaw14zhzJ76tfvdsFHtPDdkev
def get_pulls(plot_all, layers=12, cov=0.001, scatter_sigma_rad=0.01):

    ## Initialising
    n_update = 50

    # Parameter
    # start_params = np.array([-0.08146, 0.79576, -0.08657])  # [y, k, theta1]
    start_params = np.array([0.0, 0.0, 0.0])  # [y, k, theta1]
    # start_params = np.array([-4.42001868, 1.05611777, -1.19678728])  # [y, k, theta1]
    # start_params = np.array([0.0, 0.79, 0.0])  # [y, k, theta1]
    # start_params = np.array([12.345, 1, 0.0, 0.0])  # [y, k, theta1, theta2]
    delta_params = np.zeros_like(start_params)

    # Geometry
    # geo_layers = np.array(
    #     [2, 3, 5, 5.5, 5.7, 6.5, 7, 10, 12, 14, 14.01, 17, 19, 20, 21, 22, 23, 24, 25]
    # )
    geo_layers = np.random.uniform(1, layers, layers)
    geo_layers.sort()
    cov_meas = np.ones_like(geo_layers) * cov
    geo_scatter_sigma = np.zeros_like(geo_layers)
    geo_scatter_sigma[[np.random.randint(3, layers - 3)]] = scatter_sigma_rad
    # geo_scatter_sigma[[5, 9]] = scatter_sigma_rad

    # detector_layers = geo_layers[geo_scatter_sigma == 0]
    # scatter_layers = geo_layers[geo_scatter_sigma != 0]
    # print(scatter_layers)

    # Hits
    true_params = [0, 0.79]
    scatter_params = np.array([], dtype=float)
    for s_sig in geo_scatter_sigma:
        if s_sig:
            theta = c2u.scatter(s_sig)
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

    updated_params = start_params

    if plot_all:
        print(f"\ntheta = {theta}")
        print(f"measurments_all = {measurments_all}")

    chi2old = np.inf
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

            if geo_scatter_sigma[g]:  # Scatter layer
                x_s[i_s] = x

                ai = np.zeros([len(start_params), len(start_params)])
                ai[2 + i_s, 2 + i_s] = 1 / geo_scatter_sigma[g] ** 2

                bi = np.zeros_like(start_params)
                bi[2 + i_s] = -updated_params[2 + i_s] / geo_scatter_sigma[g] ** 2

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
                    logging.warning(f"i_s = {i_s} should not happen")

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

        delta_params = np.linalg.solve(a, b.transpose())

        if plot_all:
            c2u.plot_current_state(
                updated_params,
                true_params,
                a,
                np.linalg.inv(a),
                measurments_all,
                geo_layers,
                geo_scatter_sigma,
                predicted_hits,
                measurments_raw,
                n,
                "",
                "2D-Fit [y,phi]",
            )

        delta_chi2 = 1e-6
        if abs(delta_params).sum() < delta_chi2:
            logging.debug(
                f"break: 'abs(delta_params).sum() < {delta_chi2}', max updates = {n}"
            )
            break

        if chi2sum > chi2old * (1 + 1e-4):
            logging.info(f"break: 'chi2sum > chi2old', max updates = {n}")
            # break

        chi2old = chi2sum

    updated_params[1] = c2u.map_angle_to_right_half(updated_params[1], 0)
    updated_params[2] = c2u.map_angle_to_right_half(
        updated_params[2], updated_params[1]
    )

    updated_cov = np.linalg.inv(a)
    params_res, params_pulls = c2u.calc_res_pulls(
        updated_params, true_params, updated_cov
    )

    if plot_all:
        c2u.plot_current_state(
            updated_params,
            true_params,
            a,
            updated_cov,
            measurments_all,
            geo_layers,
            geo_scatter_sigma,
            predicted_hits,
            measurments_raw,
            n,
            params_pulls,
            "2D-Fit [y,phi]",
        )
        print(f"delta_params = {delta_params}")

    return params_res, params_pulls, chi2sum


logging.getLogger().setLevel(logging.INFO)
np.random.seed(10)
draws = 10000
layers = 10
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
    p_res, p_pulls, c2s = get_pulls((d < 0 or d == -1), layers)

    y_pul.append(p_pulls[0])
    k_pul.append(p_pulls[1])
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

plt.plot(t1_res, ".")
# plt.ylim(-0.5, 0.2)
plt.show()
