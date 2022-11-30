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
def straightLineCovarianceTransport(covStartParams,startParams,hPosition):
    assert(len(startParams) == 2)
    assert(covStartParams[0,1] == 0)
    assert(covStartParams[1,0] == 0)
    
    s_x = covStartParams[0,0]
    s_phi = covStartParams[1,1]
    a = hPosition/np.cos(startParams[1])
    
    transportedCov = np.array([
        [s_x+a**2*s_phi,    a*s_phi],
        [a*s_phi,           s_phi],
        ])

    return transportedCov

def projectMatrix(M,proj):
    return np.matmul(proj,np.matmul(M,proj.transpose()))




## Parameter
startParams = np.array([
    2., # x
    -0.5, # phi
    ])

covStartParams = np.array([
    [.2,0],
    [0,0.1],
    ])

nUpdate = 5
maxHorizontal = 11
maxVertical = 8
updatedParams = startParams
deltaParams = np.zeros_like(updatedParams)
detectorLayers = np.array([2,3,5,5.5,5.7,10])
measurments = np.array([1,2,3,4,3,5.5])
proj = np.array([[1,0]]) # projects onto x



## Set up plotting
fig, ax = plt.subplots()


## Iterating and updating parameters
for iUpdate in range(nUpdate):
    print(f"\nStarting iteration {iUpdate}")
    
    updatedParams = updatedParams + deltaParams
    a = np.zeros_like(covStartParams)
    b = np.zeros_like(startParams)
    chi2sum = 0
    
    # Iterate over surfaces
    for d in range(len(detectorLayers)):
        h = detectorLayers[d]
        propagatedCov = straightLineCovarianceTransport(covStartParams,updatedParams,h)
        Vi = projectMatrix(propagatedCov,proj)[0]*np.random.rand(1)*10
        ri = measurments[d] - straightLinePropagator(updatedParams,h)
        chi2i = chi2_1D(Vi,ri)
        chi2sum += chi2i
        
        print(f"Surface {d}: ri = {ri:.3f}, chi2i = {chi2i[0]:.3f}")
        
        ai = 1/Vi*np.array([
            [1,h/np.cos(updatedParams[1])**2],
            [h/np.cos(updatedParams[1])**2,(h/np.cos(updatedParams[1])**2)**2],
            ])
        bi = ri/Vi*np.array([1,h/np.cos(updatedParams[1])**2])
        
        a += ai
        b += bi
        
    print(f"chi2sum =  {chi2sum[0]:.3f}")
    
    #deltaParams = np.array([0.1,-0.3]) # WARNING experimental. add real update!

    deltaParams = np.linalg.solve(a,b.transpose())
    
    # Plot Updated Trajectory
    updatedTraj = np.array([
        [0, maxHorizontal],
        straightLinePropagator(updatedParams,[0,maxHorizontal]),
        ])    
    ax.plot(0,updatedParams[0],'xc')
    ax.plot(updatedTraj[0], updatedTraj[1],'c-')



## continue plotting
# Start Trajectory
startTraj = np.array([
    [0, maxHorizontal],
    straightLinePropagator(startParams,[0,maxHorizontal]),
    ])

ax.plot(0,startParams[0],'xr')
ax.plot(startTraj[0], startTraj[1],'r-', label='Start Trajectory')

# Detector
for d in range(len(detectorLayers)):
    ax.plot([detectorLayers[d],detectorLayers[d]],[-maxVertical,maxVertical],'g-')
    ax.plot(detectorLayers[d],measurments[d],'gx')
    
# Updated Trajectory
updatedTraj = np.array([
    [0, maxHorizontal],
    straightLinePropagator(updatedParams,[0,maxHorizontal]),
    ])

ax.plot(0,updatedParams[0],'xb')
ax.plot(updatedTraj[0], updatedTraj[1],'b-', label='Final Trajectory')

ax.set(xlabel='horizontal', ylabel='x',
       title='2D-Fit')
ax.legend()

#fig.savefig("test.png")
plt.show()
