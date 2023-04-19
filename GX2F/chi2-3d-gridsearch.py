# Alexander J. Pfleger
# 2023-04-13
#
# Gridsearch

import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib import ticker, cm
# from matplotlib.ticker import LinearLocator

import chi2_utilities as c2u
import propagators

# Initialising
cov=0.1
scatter_sigma_rad=0.0316

# Geometry
geo_layers = np.array(
    [2, 3, 5, 5.5, 5.7, 6.5, 7, 10, 12, 14, 14.01, 17, 19, 20, 21, 22, 23, 24, 25]
)
cov_meas = np.ones_like(geo_layers) * cov
geo_scatter_sigma = np.zeros_like(geo_layers)
geo_scatter_sigma[[10]] = scatter_sigma_rad

# measurments_raw = np.array([ 2.01849258,  3.0277388,   5.04623144,  5.55085459,
#                              5.75270384,  6.56010087,  7.06472402, 10.09246288,
#                             12.11095546, 14.12944804, 14.1395405,  16.72095812,
#                             18.44765888, 19.31100925, 20.17435963, 21.03771,
#                             21.90106038, 22.76441075, 23.62776113])

measurments_raw = np.array([ 1.2559942,   2.71279278,  5.52939274,  5.83925736,
                             5.85874512,  6.43183048,  7.02965994, 10.26960115,
                            11.99163701, 14.5388773,  14.61684614, 17.0698621,
                            18.26044612, 19.26869784, 20.15545842, 21.13038446,
                            21.45491631, 23.16104188, 23.34623096])



# alpha = np.array([0, 0, 0], dtype=float)
# delta_alpha = np.array([100, np.pi/2*0.9, np.pi/2*0.9], dtype=float)

alpha = np.array([-4.42033707,  1.05613163, -1.19687986], dtype=float)
delta_alpha = np.array([0.1, 0.01, 0.01], dtype=float)

steps = 10
alpha_old = alpha + np.inf

print(alpha + delta_alpha)
while sum(abs(alpha_old / alpha - 1)) > 1e-7:
    # chi2min = np.inf
    chi2min = 0
    alpha_old = alpha.copy()
    alpha_vec = np.zeros([3,steps], dtype=float)
    for a_dim in range(3):
        alpha_vec[a_dim] = np.linspace(alpha[a_dim]-delta_alpha[a_dim], alpha[a_dim]+delta_alpha[a_dim], steps)
    print(alpha_vec.T[-1])
        
    for y in alpha_vec[0]:
        for phi in alpha_vec[1]:
            for theta in alpha_vec[2]:
    
                if abs(theta + phi) >= np.pi/2:
                    continue
    
                chi2sum = 0
        
                fit_params = np.array([y, phi, theta])
                
                predicted_hits = propagators.straight_line_propagator_stepwise_2D_scatter(
                    fit_params, geo_layers, geo_scatter_sigma, "phi"
                )
                
                residuals = measurments_raw - predicted_hits
                
                chi2sum = c2u.chi2_1D(cov, residuals).sum()
                chi2sum += c2u.chi2_1D(scatter_sigma_rad ** 2, fit_params[2])
                
                if chi2sum > chi2min:
                    chi2min = chi2sum
                    alpha = fit_params
                    
                    # print(f"chi2 = {chi2sum:0.5f} @ {y:0.5f} | {phi:0.5f} | {theta:0.5f}")
                    
    delta_alpha = delta_alpha / np.sqrt(steps)
    
    print(alpha)
