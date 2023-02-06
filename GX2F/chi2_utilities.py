import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import chi2


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


def add_traj_to_plot(
    ax, params, maxHorizontal, propagator, color="b", label_text="", style="-"
):
    traj = np.array(
        [
            [0, maxHorizontal],
            propagator(params, [0, maxHorizontal]),
        ]
    )
    try:
        ax.plot(0, params[0], "x" + color)
    except:
        ax.plot(0, params, "x" + color)
    ax.plot(traj[0], traj[1], color + style, label=label_text)


def plot_pull_distribution(pulls, title):
    bins = int(np.sqrt(len(pulls)))
    mu, std = norm.fit(pulls)

    fig, ax = plt.subplots()
    plt.hist(pulls, bins=bins, density=True)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 201)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, "k")
    plt.title(f"{title}: $\mu$ = {mu:.3f}, $\sigma$ = {std:.3f}")
    plt.show()


def plot_chi2_distribution(chi2sum, title):
    bins = int(np.sqrt(len(chi2sum)))
    df, loc, scale = chi2.fit(chi2sum)

    fig, ax = plt.subplots()
    plt.hist(chi2sum, bins=bins, density=True)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 201)
    p = chi2.pdf(x, df, loc, scale)
    plt.plot(x, p, "k")
    plt.title(f"{title}: k = {df:.3f}, loc = {loc:.3f}, scale = {scale:.3f}")
    plt.show()