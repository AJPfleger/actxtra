import numpy as np


def generate_hits(geometry, true_params, propagator, cov=0.1, smearing=True):
    measurments_raw = propagator(true_params, geometry)

    cov_meas = [cov] * len(measurments_raw)

    measurments = []
    for mi in range(len(measurments_raw)):
        if smearing:
            m = np.random.normal(measurments_raw[mi], np.sqrt(cov_meas[mi]))
        else:
            m = measurments_raw[mi]
        measurments.append(m)

    return measurments, cov_meas, measurments_raw


def chi2_1D(V, r):
    return r ** 2 / V  # r * (1 / V) * r
