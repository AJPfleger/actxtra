# Alexander J. Pfleger
# 2023-03-20
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 2D case
# params are [y, phi, theta1, theta2]
# scattering is implemented via active surfaces, that give a gaussian smearing
# to the direction. This example includes 2 scattering surfaces

import numpy as np
import logging

import chi2_utilities as c2u
import propagators


def get_chi2sum(
    params, measurments_raw, geo_layers, geo_scatter_sigma, cov, scatter_sigma_rad
):

    chi2sum = 0

    predicted_hits = propagators.straight_line_propagator_stepwise_2D_scatter(
        params, geo_layers, geo_scatter_sigma, "phi"
    )

    residuals = measurments_raw - predicted_hits
    residuals = residuals[geo_scatter_sigma == 0]

    chi2sum = c2u.chi2_1D(cov, residuals).sum()
    chi2sum += c2u.chi2_1D(scatter_sigma_rad ** 2, params[2])
    chi2sum += c2u.chi2_1D(scatter_sigma_rad ** 2, params[3])

    return chi2sum


def gridsearch(
    measurments_raw, geo_layers, geo_scatter_sigma, cov_meas, cov, scatter_sigma_rad
):

    alpha = np.array([0, 0, 0, 0], dtype=float)
    delta_alpha = np.array(
        [10, np.pi / 2 * 0.9, np.pi / 2 * 0.9, np.pi / 2 * 0.9], dtype=float
    )

    steps = 7
    alpha_old = alpha + np.inf
    it_count = 0

    while (
        sum(abs(alpha_old / (alpha + 1e-15) - 1)) > 1e-4
    ):  # add 1e-15 to prevent divison by zero
        it_count += 1
        logging.debug(f"it_count: {it_count}")
        chi2min = np.inf
        alpha_old = alpha.copy()
        alpha_vec = np.zeros([4, steps], dtype=float)
        for a_dim in range(4):
            alpha_vec[a_dim] = np.linspace(
                alpha[a_dim] - delta_alpha[a_dim],
                alpha[a_dim] + delta_alpha[a_dim],
                steps,
            )

        for y in alpha_vec[0]:
            for phi in alpha_vec[1]:
                for theta1 in alpha_vec[2]:
                    for theta2 in alpha_vec[3]:

                        if abs(phi + theta1 + theta2) >= np.pi / 2:
                            continue

                        fit_params = np.array([y, phi, theta1, theta2])
                        chi2sum = get_chi2sum(
                            fit_params,
                            measurments_raw,
                            geo_layers,
                            geo_scatter_sigma,
                            cov,
                            scatter_sigma_rad,
                        )

                        if chi2sum < chi2min:
                            chi2min = chi2sum
                            alpha = fit_params

                            logging.debug(
                                f"chi2 = {chi2sum:0.5f} @ {y:0.5f} | {phi:0.5f} | {theta1:0.5f} | {theta2:0.5f}"
                            )

        delta_alpha = delta_alpha / 2

    logging.info(f"gridsearch finished after {it_count} iterations")

    return alpha, chi2min


# scatter-angle for 1 GeV pion with 300 um silicon is about 0.8 mrad
# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwizosSqoLv-AhXUGogKHQYBBKkQFnoECCEQAQ&url=https%3A%2F%2Fe-publishing.cern.ch%2Findex.php%2FCYRSP%2Farticle%2Fdownload%2F534%2F396%2F1738&usg=AOvVaw14zhzJ76tfvdsFHtPDdkev
def get_pulls(plot_all, layers=12, cov=0.001, scatter_sigma_rad=0.05):

    ## Initialising
    n_update = 50

    # Parameter
    start_params = np.array([0.0, 0.0, 0.0, 0.0])  # [y, k, theta1, theta2]
    delta_params = np.zeros_like(start_params)

    # Geometry
    geo_layers = np.array(
        [
            2,
            3,
            5,
            5.5,
            5.7,
            6.5,
            7,
            10,
            12,
            14,
            14.01,
            17,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            26.3,
            26.7,
            27,
            30,
            32.5,
            33,
            33.5,
            34,
            37,
            40,
        ]
    )
    # geo_layers = np.random.uniform(1, layers, layers)
    geo_layers.sort()
    cov_meas = np.ones_like(geo_layers) * cov
    geo_scatter_sigma = np.zeros_like(geo_layers)
    # geo_scatter_sigma[[np.random.randint(3, layers - 3, 2)]] = scatter_sigma_rad
    geo_scatter_sigma[[11, 20]] = scatter_sigma_rad

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
        updated_params[3] = c2u.map_angle_to_right_half(
            updated_params[3], updated_params[1] + updated_params[2]
        )

        predicted_hits = propagators.straight_line_propagator_stepwise_2D_scatter_yphi(
            updated_params, geo_layers, geo_scatter_sigma
        )

        # Iterate over surfaces
        a = np.zeros([len(start_params), len(start_params)])
        b = np.zeros_like(start_params)
        chi2sum = 0
        i_s = 0

        phi = updated_params[1]
        theta1 = updated_params[2]
        theta2 = updated_params[3]
        cosphi2 = np.cos(phi) ** 2
        cosphitheta12 = np.cos(phi + theta1) ** 2
        cosphitheta22 = np.cos(phi + theta1 + theta2) ** 2

        x_s = np.array([0, 0])

        residuals = measurments_all - predicted_hits

        for g in range(len(geo_layers)):
            x = geo_layers[g]

            if geo_scatter_sigma[g]:  # Scatter layer
                x_s[i_s] = x

                ai = np.zeros([len(start_params), len(start_params)])
                ai[2 + i_s, 2 + i_s] = 1 / geo_scatter_sigma[g] ** 2

                bi = np.zeros_like(start_params)
                bi[2 + i_s] = -updated_params[2 + i_s] / geo_scatter_sigma[g] ** 2

                chi2sum += c2u.chi2_1D(
                    geo_scatter_sigma[g] ** 2, updated_params[2 + i_s]
                )

                i_s += 1
            else:  # Detector layer
                Vi = cov_meas[g]
                ri = residuals[g]
                chi2sum += c2u.chi2_1D(Vi, ri)

                if i_s == 0:
                    dydy0 = 1
                    dydp0 = x / cosphi2
                    dydt1 = 0
                    dydt2 = 0
                elif i_s == 1:
                    dydy0 = 1
                    dydp0 = x_s[0] / cosphi2 + (x - x_s[0]) / cosphitheta12
                    dydt1 = (x - x_s[0]) / cosphitheta12
                    dydt2 = 0
                elif i_s == 2:
                    dydy0 = 1
                    dydp0 = (
                        x_s[0] / cosphi2
                        + (x_s[1] - x_s[0]) / cosphitheta12
                        + (x - x_s[1]) / cosphitheta22
                    )
                    dydt1 = (x_s[1] - x_s[0]) / cosphitheta12 + (
                        x - x_s[1]
                    ) / cosphitheta22
                    dydt2 = (x - x_s[1]) / cosphitheta22
                else:
                    logging.warning(f"i_s = {i_s} should not happen")

                abi_vec = np.array([[dydy0, dydp0, dydt1, dydt2]])

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

        if chi2sum > chi2old * (1 + 1e-4) and n != 1:
            logging.info(f"'chi2sum > chi2old', max updates = {n}")
            logging.info(
                f"chi2sum = {chi2sum}\n"
                f"chi2old = {chi2old}\n"
                f"rel. error = {chi2sum / chi2old - 1}\n"
            )
            # break
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

        chi2old = chi2sum
        # print(f"chi2sum --- {chi2sum}")

    updated_params[1] = c2u.map_angle_to_right_half(updated_params[1], 0)
    updated_params[2] = c2u.map_angle_to_right_half(
        updated_params[2], updated_params[1]
    )
    updated_params[3] = c2u.map_angle_to_right_half(
        updated_params[3], updated_params[1] + updated_params[2]
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

        params_grid, chi2_grid = gridsearch(
            measurments_all,
            geo_layers,
            geo_scatter_sigma,
            cov_meas,
            cov,
            scatter_sigma_rad,
        )
        print(f"params_grid   = {params_grid}\tchi2: {chi2_grid:0.3f}")
        print(f"params_normal = {updated_params}\tchi2: {chi2sum:0.3f}\n")

    return params_res, params_pulls, chi2sum


logging.getLogger().setLevel(logging.INFO)
np.random.seed(100)
draws = 10000
layers = 30
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
    logging.debug(f"\ndraw {d}")
    p_res, p_pulls, c2s = get_pulls((d < -1 or d == -1), layers)

    y_pul.append(p_pulls[0])
    k_pul.append(p_pulls[1])
    t1_pul.append(p_pulls[2])
    t2_pul.append(p_pulls[3])
    y_res.append(p_res[0])
    k_res.append(p_res[1])
    t1_res.append(p_res[2])
    t2_res.append(p_res[3])
    chi2sum.append(c2s)

c2u.plot_pull_distribution(y_res, f"y_res ({layers} hits)")
c2u.plot_pull_distribution(k_res, f"phi_res({layers} hits)")
c2u.plot_pull_distribution(t1_res, f"theta1_res ({layers} hits)")
c2u.plot_pull_distribution(t2_res, f"theta2_res ({layers} hits)")

c2u.plot_pull_distribution(y_pul, f"y_pulls ({layers} hits)")
c2u.plot_pull_distribution(k_pul, f"phi_pulls ({layers} hits)")
c2u.plot_pull_distribution(t1_pul, f"theta1_pulls ({layers} hits)")
c2u.plot_pull_distribution(t2_pul, f"theta2_pulls ({layers} hits)")

c2u.plot_chi2_distribution(chi2sum, f"$\chi^2$ ([y,k], {layers} hits)")
