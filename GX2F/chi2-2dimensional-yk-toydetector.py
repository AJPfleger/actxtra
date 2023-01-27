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
from scipy.stats import norm

import chi2_utilities as c2u


#checked
def straight_line_propagator(params, x_vec):
    y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

    return y_vec


def get_pulls(plot_all, layers=12, cov=0.1):
    
    ## Initialising
    n_update = 15

    # Parameter
    start_params = np.array([0.0, 3])  # [y, k]
    delta_params = np.zeros_like(start_params)
    # updatedCov = np.array([
    #     [3,0],
    #     [0,0.1],
    #     ])

    # Geometry
    # detectorLayers = np.array([2, 3, 5, 5.5, 5.7, 6.5])  # ,10,50,98,99,100])
    detectorLayers = np.random.uniform(1, 12, layers)
    
    # Hits
    true_params = [12.345, 1]
    # true_params = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments, cov_meas, measurments_raw = c2u.generate_hits(
        detectorLayers, true_params, straight_line_propagator, cov, True
    )
    # proj = np.array([[1,0]]) # projects onto x

    # start_params = true_params
    updated_params = start_params

    if plot_all:
        maxHorizontal = max(detectorLayers) + 1
        maxVertical = 40
        ## Set up plotting
        fig, ax = plt.subplots()

        def add_traj_to_plot(
            params, color="b", label_text="", style="-", maxHorizontal=maxHorizontal
        ):
            traj = np.array(
                [
                    [0, maxHorizontal],
                    straight_line_propagator(params, [0, maxHorizontal]),
                ]
            )

            ax.plot(0, params[0], "x" + color)
            ax.plot(traj[0], traj[1], color + style, label=label_text)
    
    ## Iterating and updating parameters
    for _ in range(n_update):

        updated_params = updated_params + delta_params
        
        # Iterate over surfaces
        a = np.zeros([2, 2])
        b = np.zeros_like(start_params)
        chi2sum = 0
        for d in range(len(detectorLayers)):
            h = detectorLayers[d]
            # propagatedCov = straightLineCovarianceTransport(updatedCov,updated_params,h)
            Vi = cov_meas[d]
            ri = measurments[d] - straight_line_propagator(updated_params, h)
            chi2sum += c2u.chi2_1D(Vi, ri)
            
            ai = (
                1 / Vi * np.array([
                        [1,      h],
                        [h, h**2],
                        ])
            )
            bi = ri / Vi * np.array([1, h])

            a += ai
            b += bi

        delta_params = np.linalg.solve(a, b.transpose())

    updatedCov = np.linalg.inv(a)

    y_res, k_res = updated_params - true_params
    y_pull = y_res / np.sqrt(updatedCov[0][0])
    k_pull = k_res / np.sqrt(updatedCov[1][1])

    if plot_all:
        print(f"updated_params: {updated_params}")
        print(f"true_params: {true_params}")
        print(f"diff: {updated_params - true_params}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updatedCov:\n{updatedCov}")
        print(f"pulls: {y_pull}, {k_pull}")
        print("\n")

        # continue plotting
        # Detector
        for d in range(len(detectorLayers)):
            ax.plot(
                [detectorLayers[d], detectorLayers[d]],
                [-maxVertical, maxVertical],
                "g-",
            )
            ax.plot(detectorLayers[d], measurments[d], "gx")

        # Updated Trajectory

        # Plot Trajectoris
        add_traj_to_plot(start_params, "r", "Start Trajectory", "-")
        add_traj_to_plot(updated_params, "b", "Final Trajectory", "-")
        add_traj_to_plot(true_params, "k", "Unsmeared True Trajectory", "-.")

        ax.set(xlabel="horizontal", ylabel="x", title="2D-Fit")
        ax.legend()

        # fig.savefig("test.png")
        plt.show()

    return y_pull, k_pull


draws = 10000
bins = int(np.sqrt(draws))
y_pulls = []
k_pulls = []
for d in range(draws):
    y_p, k_p = get_pulls(d < 5)
    y_pulls.append(y_p)
    k_pulls.append(k_p)




mu, std = norm.fit(y_pulls)
fig, ax = plt.subplots()
ax.hist(y_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
ax.plot(x, p, "k")
ax.set(title=f"y_pulls: mu = {mu:.5f}, std = {std:.5f}")
fig.savefig("yk-y_pulls.png")
plt.show()



mu, std = norm.fit(k_pulls)
fig, ax = plt.subplots()
ax.hist(k_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
ax.plot(x, p, "k")
ax.set(title=f"k_pulls: mu = {mu:.5f}, std = {std:.5f}")
fig.savefig("yk-k_pulls.png")
plt.show()

