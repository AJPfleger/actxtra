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


def chi2_1D(V,r):
    return r*(1/V)*r


# def derive1Chi2(H,V,r):
#     return -2*H.transpose()*np.linalg.inv(V)*r

# def derive2Chi2(H,V):
#     return 2*H.transpose()*np.linalg.inv(V)*H


def straightLinePropagator(startParams,hPosition):
    xVec = np.ones_like(hPosition)*startParams[0]
    phiVec = np.ones_like(hPosition)*startParams[1]
    return xVec+hPosition*np.tan(phiVec)


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


def projectMatrix(M,proj):
    return np.matmul(proj,np.matmul(M,proj.transpose()))


def generateHits(geometry, start_params_hits):
    
    measurments_raw = straightLinePropagator(start_params_hits, geometry)
    #measurments = np.array([1,2,3,4,3,5.5])
    #measurments = measurments_raw
    
    cov_meas = [ 0.1 for _ in range(len(measurments_raw))]

    measurments = []
    for mi in range(len(measurments_raw)):
        m = np.random.normal(measurments_raw[mi], cov_meas[mi])
        measurments.append(m)

    return measurments, cov_meas, measurments_raw


# def add_traj_to_plot(params, color="b", label_text="", style="-", maxHorizontal=11):
#     traj = np.array([
#         [0, maxHorizontal],
#         straightLinePropagator(params,[0,maxHorizontal]),
#         ])
    
#     ax.plot(0, params[0],"x"+color)
#     ax.plot(traj[0], traj[1], color+style, label=label_text)
    
    
def getPulls(plot_all):

    ## Parameter
    startParams = np.array([
        0., # x
        0., # phi
        ])
    
    
    # covStartParams = np.array([
    #     [3,0],
    #     [0,0.1],
    #     ])

    nUpdate = 5
    updatedParams = startParams
    deltaParams = np.zeros_like(updatedParams)
    detectorLayers = np.array([2,3,5,5.5,5.7,6.5])#,10,50,98,99,100])
    start_params_hits = [5.1, -0.9]
    measurments, cov_meas, measurments_raw = generateHits(detectorLayers, start_params_hits)
    # proj = np.array([[1,0]]) # projects onto x
    
    startParams = start_params_hits
    
    if plot_all:
        maxHorizontal = max(detectorLayers) + 1
        maxVertical = 40
        ## Set up plotting
        fig, ax = plt.subplots()
        def add_traj_to_plot(params, color="b", label_text="", style="-", maxHorizontal=maxHorizontal):
            traj = np.array([
                [0, maxHorizontal],
                straightLinePropagator(params,[0,maxHorizontal]),
                ])
            
            ax.plot(0, params[0],"x"+color)
            ax.plot(traj[0], traj[1], color+style, label=label_text)
    
    # updatedCov = covStartParams
    ## Iterating and updating parameters
    for iUpdate in range(nUpdate):
        # print(f"\nStarting iteration {iUpdate}")
        
        updatedParams = updatedParams + deltaParams
        a = np.zeros([2,2])
        b = np.zeros_like(startParams)
        chi2sum = 0
        
        # Iterate over surfaces
        for d in range(len(detectorLayers)):
            h = detectorLayers[d]
            # propagatedCov = straightLineCovarianceTransport(updatedCov,updatedParams,h)
            # Vi = projectMatrix(propagatedCov,proj)[0]
            Vi = cov_meas[d]
            ri = measurments[d] - straightLinePropagator(updatedParams,h)
            chi2i = chi2_1D(Vi,ri)
            chi2sum += chi2i
            
            # print(f"Surface {d}: ri = {ri:.3f}, chi2i = {chi2i:.3f}")
            
            ai = 1/Vi*np.array([
                [1, h/np.cos(updatedParams[1])**2],
                [h/np.cos(updatedParams[1])**2, (h/np.cos(updatedParams[1])**2)**2],
                ])
            bi = ri/Vi*np.array([1, h/np.cos(updatedParams[1])**2])
            
            a += ai
            b += bi
        # updatedCov = np.linalg.inv(a)
        # print(f"chi2sum =  {chi2sum:.3f}")
        # print(f'parameter: {updatedParams}')
        
        #deltaParams = np.array([0.1,-0.3]) # WARNING experimental. add real update!
    
        deltaParams = np.linalg.solve(a,b.transpose())
        
        # Plot Updated Trajectory
        # add_traj_to_plot(updatedParams, "c", "", "-")
    
    # print(f'a:\n{a}')
    updatedCov = np.linalg.inv(a)
    # print(f'updatedCov:\n{updatedCov}')
    
    y_pull, phi_pull = (updatedParams - start_params_hits) #/ [updatedCov[0][0], updatedCov[1][1]]
    
    
    if plot_all:
        print(f'updatedParams: {updatedParams}')
        print(f'start_params_hits: {start_params_hits}')
        print(f'diff: {updatedParams - start_params_hits}')
        print(f'updatedCov:\n{updatedCov}')
        print(f'pulls: {y_pull}, {phi_pull}')
        print("\n")
        
        # continue plotting
        # Detector
        for d in range(len(detectorLayers)):
            ax.plot([detectorLayers[d],detectorLayers[d]],[-maxVertical,maxVertical],'g-')
            ax.plot(detectorLayers[d],measurments[d],'gx')
            
        # Updated Trajectory
        
        # Plot Trajectoris
        add_traj_to_plot(startParams, "r", "Start Trajectory", "-")
        add_traj_to_plot(updatedParams, "b", "Final Trajectory", "-")
        add_traj_to_plot(start_params_hits, "k", "Unsmeared True Trajectory", "-.")
        
        ax.set(xlabel='horizontal', ylabel='x',
                title='2D-Fit')
        ax.legend()
        
        #fig.savefig("test.png")
        plt.show()
    
    return y_pull, phi_pull

draws = 1000
bins = int(np.sqrt(draws))
y_pulls = []
phi_pulls = []
for d in range(draws):
    y_p, phi_p = getPulls(d<5)
    y_pulls.append(y_p)
    phi_pulls.append(phi_p)

from scipy.stats import norm

mu, std = norm.fit(y_pulls)
plt.hist(y_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k')
plt.title(f"y_pulls: mu = {mu:.5f}, std = {std:.5f}")
plt.show()

mu, std = norm.fit(phi_pulls)
plt.hist(phi_pulls, bins=bins, density=True)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k')
plt.title(f"phi_pulls: mu = {mu:.5f}, std = {std:.5f}")
plt.show()
