# Alexander J. Pfleger
# 2022-11-25
#
# Example to test the functionality of the chi2-algorithm
# propagation is done by a straight line
# surfaces are infinitly large parallel planes
# 1D case

import numpy as np
import matplotlib.pyplot as plt


# def residual(m,H,x):
#     return m - H*x


# def chi2(V,r):
#     return r.transpose()*np.linalg.inv(V)*r


def chi2_1D(V, r):
    return r * (1 / V) * r


# def derive1Chi2(H,V,r):
#     return -2*H.transpose()*np.linalg.inv(V)*r


# def derive2Chi2(H,V):
#     return 2*H.transpose()*np.linalg.inv(V)*H

#checked
def straight_line_propagator(params, x_vec):

    if not (-np.pi/2 < params[1] < np.pi/2):
        print(f"ERROR straight_line_propagator: phi {params[1]} is out of bounds")
    
    y_vec = np.ones_like(x_vec) * params[0] + x_vec * np.ones_like(x_vec) * params[1]

    return y_vec


# works only for x-phi-case with diagonal covariance
# def straightLineCovarianceTransport(covStartParams,startParams,hPosition):
#     assert(len(startParams) == 2)
#     #assert(covStartParams[0,1] == 0)
#     #assert(covStartParams[1,0] == 0)

#     s_x = covStartParams[0,0]
#     s_phi = covStartParams[1,1]
#     a = hPosition/np.cos(startParams[1])

#     transportedCov = np.array([
#         [s_x+a**2*s_phi,    a*s_phi],
#         [a*s_phi,           s_phi],
#         ])

#     return transportedCov


def projectMatrix(M, proj):
    return np.matmul(proj, np.matmul(M, proj.transpose()))

#checked
def generate_hits(geometry, start_params_hits):

    measurments_raw = straight_line_propagator(start_params_hits, geometry)

    cov_meas = [0.1] * len(measurments_raw)

    measurments = []
    for mi in range(len(measurments_raw)):
        m = np.random.normal(measurments_raw[mi], np.sqrt(cov_meas[mi]))
        # m = measurments_raw[mi]
        measurments.append(m)

    return measurments, cov_meas, measurments_raw


def getPulls(plot_all):
    
    ## Initialising
    n_update = 15

    # Parameter
    start_params = np.array([0.0, 0.0,])  # [y, phi]
    delta_params = np.zeros_like(start_params)
    # updatedCov = np.array([
    #     [3,0],
    #     [0,0.1],
    #     ])

    # Geometry
    detectorLayers = np.array([2, 3, 5])  #, 5.5, 5.7, 6.5])  # ,10,50,98,99,100])
    
    # Hits
    start_params_hits = [5.1, -0.9]
    # start_params_hits = [np.random.uniform(-9,9), np.random.uniform(-0.9,0.9)]
    measurments, cov_meas, measurments_raw = generate_hits(
        detectorLayers, start_params_hits
    )
    # proj = np.array([[1,0]]) # projects onto x

    # start_params = start_params_hits
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
            # Vi = projectMatrix(updatedCov,proj)[0]
            Vi = cov_meas[d]
            ri = measurments[d] - straight_line_propagator(updated_params, h)
            chi2sum += chi2_1D(Vi, ri)
            
            h_cos2 = h / np.cos(updated_params[1]) ** 2
            d2chi_dphi2 = h_cos2 ** 2 #- ri*h_cos2*2*np.tan(updated_params[1])
            ai = (
                1 / Vi * np.array([
                        [1,      h_cos2],
                        [h_cos2, d2chi_dphi2],
                        ])
            )
            bi = ri / Vi * np.array([1, h_cos2])

            a += ai
            b += bi
        # updatedCov = np.linalg.inv(a)
        # print(f"chi2sum =  {chi2sum:.3f}")
        # print(f'parameter: {updated_params}')

        # delta_params = np.array([0.1,-0.3]) # WARNING experimental. add real update!

        delta_params = np.linalg.solve(a, b.transpose())

        # Plot Updated Trajectory
        # add_traj_to_plot(updated_params, "c", "", "-")
    
    # print(f'a:\n{a}')
    updatedCov = np.linalg.inv(a)
    # print(f'updatedCov:\n{updatedCov}')

    y_pull, phi_pull = (
        updated_params - start_params_hits
    )  / [np.sqrt(updatedCov[0][0]), np.sqrt(updatedCov[1][1])]

    if plot_all:
        print(f"updated_params: {updated_params}")
        print(f"start_params_hits: {start_params_hits}")
        print(f"diff: {updated_params - start_params_hits}")
        print(f"a:\n{a}")
        print(f"cov_meas: {cov_meas}")
        print(f"updatedCov:\n{updatedCov}")
        print(f"pulls: {y_pull}, {phi_pull}")
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
        add_traj_to_plot(start_params_hits, "k", "Unsmeared True Trajectory", "-.")

        ax.set(xlabel="horizontal", ylabel="x", title="2D-Fit")
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
