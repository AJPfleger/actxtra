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


def fit_func(x, y0, phi0):
    return y0 + np.tan(phi0) * x


#checked
def straight_line_propagator(params, x_vec):

    if not (-np.pi/2 < params[1] < np.pi/2):
        print(f"ERROR straight_line_propagator: phi {params[1]} is out of bounds")
        # if -np.pi/2 < params[1]:
        #     params[1] = -np.pi/2*0.999
        # else:
        #     params[1] = np.pi/2*0.999
        while not (-np.pi/2 < params[1] < np.pi/2):
            if params[1] < 0:
                params[1] += np.pi
            else:
                params[1] -= np.pi
    
    y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * np.tan(params[1])

    return y_vec


def get_pulls(plot_all, layers=5, cov=0.1, phi_true = -1):
    
    ## Initialising
    n_update = 5

    # Parameter
    start_params = np.array([0.0, 0])  # [y, phi]
    delta_params = np.zeros_like(start_params)
    # updated_cov = np.array([
    #     [3,0],
    #     [0,0.1],
    #     ])

    # Geometry
    # detector_layers = np.array([2, 3, 5, 5.5, 5.7, 6.5, 10, 50, 98, 99, 100])
    detector_layers = np.random.uniform(2, 10, layers)
    
    # Hits
    true_params = [12.345, phi_true]
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detector_layers, true_params, straight_line_propagator, 0.1, True
    )
    # proj = np.array([[1,0]]) # projects onto x
    start_params = true_params
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
            
            h_cos2 = h / np.cos(updated_params[1]) ** 2
            d2chi_dphi2 = h_cos2 ** 2 #- 2 * ri / h_cos2 ** 2 * np.tan(updated_params[1])
            ai = (
                1 / Vi * np.array([
                        [1,      h_cos2],
                        [h_cos2, d2chi_dphi2],
                        ])
            )
            bi = ri / Vi * np.array([1, h_cos2])

            a += ai
            b += bi
        # updated_cov = np.linalg.inv(a)
        # print(f"chi2sum =  {chi2sum:.3f}")
        # print(f'parameter: {updated_params}')

        # delta_params = np.array([0.1,-0.3]) # WARNING experimental. add real update!

        delta_params = np.linalg.solve(a, b.transpose())

        # Plot Updated Trajectory
        # add_traj_to_plot(updated_params, "c", "", "-")
    
    # print(f'a:\n{a}')
    updated_cov = np.linalg.inv(a)
    # print(f'updated_cov:\n{updated_cov}')

    y_res, phi_res = updated_params - true_params
    y_pull = y_res / np.sqrt(updated_cov[0][0])
    phi_pull = phi_res / np.sqrt(updated_cov[1][1])

    if plot_all:        
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updated_cov:\n{updated_cov}")
        print(f"pulls: {y_pull}, {phi_pull}")
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

    return y_pull, phi_pull


draws = 10000
bins = int(np.sqrt(draws))
y_pulls = []
phi_pulls = []
for d in range(draws):
    y_p, phi_p = getPulls(d < 5)
    y_pulls.append(y_p)
    phi_pulls.append(phi_p)

from scipy.stats import norm

mu, std = norm.fit(y_pulls)
plt.hist(y_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, "k")
plt.title(f"y_pulls: mu = {mu:.5f}, std = {std:.5f}")
plt.show()

mu, std = norm.fit(phi_pulls)
plt.hist(phi_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, "k")
plt.title(f"phi_pulls: mu = {mu:.5f}, std = {std:.5f}")
plt.show()
