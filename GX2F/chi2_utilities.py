import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import chi2
import logging


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


def scatter(sigma):
    return np.random.normal(0, sigma)


def generate_hits_scatter(
    geometry, geo_scatter_sigma, true_params, propagator, cov=0.1
):

    measurments_raw = propagator(true_params, geometry, geo_scatter_sigma)
    measurments_smeared = np.random.normal(measurments_raw, np.sqrt(cov))

    return measurments_smeared, measurments_raw


def chi2_1D(V, r):
    return r ** 2 / V  # r * (1 / V) * r


def add_traj_to_plot(
    ax, params, maxHorizontal, propagator, color="b", label_text="", style="-"
):
    traj = np.array([[0, maxHorizontal], propagator(params, [0, maxHorizontal]),])
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
    # plt.xlim(0, 100)
    # plt.savefig("second_derivatives_false.pdf")
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


def map_angle_to_right_half(beta, offset=0):
    while not (-np.pi / 2 < beta + offset < np.pi / 2):
        logging.info("angle adjustment needs to be done")
        beta += np.pi if beta + offset < 0 else -np.pi

    return beta


# only generates plot, if a proper title is set
# only saves plot, if a filename is provided
def plot_current_state(
    updated_params,
    true_params,
    a,
    updated_cov,
    measurments_all,
    geo_layers,
    geo_scatter_sigma,
    predicted_hits,
    measurments_raw,
    n="",
    params_pulls="",
    plot_title="",
    plot_filename="",
):
    if n != "":
        print(f"\nmax updates = {n}")
    print(
        f"updated_params: {updated_params}\n"
        f"true_params: {true_params}\n"
        f"diff: {updated_params - true_params}\n"
        f"a:\n{a}\n"
        f"updated_cov:\n{updated_cov}"
    )
    if n != "":
        print(f"pulls: {params_pulls}\n")
    print("\n")

    if plot_title != "":
        delta_measurments = abs(max(measurments_all) - min(measurments_all))
        min_vertical = min(measurments_all) - 0.3 * delta_measurments
        max_vertical = max(measurments_all) + 0.3 * delta_measurments

        fig, ax = plt.subplots()

        # Detector
        for d in range(len(geo_layers)):
            ax.plot(
                [geo_layers[d], geo_layers[d]],
                [min_vertical, max_vertical],
                "g-" if geo_scatter_sigma[d] == 0 else "r:",
            )

        ax.plot(geo_layers, measurments_all, "gx")

        # Trajectories
        ax.plot(
            np.append(0, geo_layers),
            np.append(updated_params[0], predicted_hits),
            "b-",
            label="Fitted Trajectory",
        )
        ax.plot(
            np.append(0, geo_layers),
            np.append(true_params[0], measurments_raw),
            "k-",
            label="Unsmeared True Trajectory",
        )

        ax.set(xlabel="x", ylabel="y", title=plot_title)
        ax.legend()

        if plot_filename != "":
            fig.savefig(plot_filename)
            logging.info(f"Plot saved as {plot_filename}")
        plt.show()


def calc_res_pulls(updated_params, true_params, params_cov):
    params_res = updated_params - true_params

    params_pulls = np.zeros_like(params_res)
    for p in range(len(params_res)):
        params_pulls[p] = params_res[p] / np.sqrt(params_cov[p][p])

    return params_res, params_pulls


def root_fit(detector_layers, measurments, cov, fit_func, true_params, angle_like=[]):
    import ROOT

    x_root = np.array(detector_layers)
    y_root = np.array(measurments)
    ex_root = x_root * 0
    ey_root = ex_root + np.sqrt(cov)
    g_root = ROOT.TGraphErrors(len(x_root), x_root, y_root, ex_root, ey_root)
    f_root = ROOT.TF1("f_root", fit_func)
    t_root = g_root.Fit(f_root, "S")

    params_res = np.zeros_like(true_params)
    params_pulls = np.zeros_like(true_params)
    for p in range(len(true_params)):
        param_p = (
            t_root.Parameter(p)
            if p not in angle_like
            else map_angle_to_right_half(t_root.Parameter(p), 0)
        )
        params_res[p] = param_p - true_params[p]
        params_pulls[p] = params_res[p] / t_root.Error(p)

    return params_res, params_pulls
