import numpy as np


def straight_line_propagator_2D_yphi(params, x_vec):
    return straight_line_propagator_2D(params, x_vec, "phi")


def straight_line_propagator_2D(params, x_vec, mode="phi"):
    y0 = np.ones_like(x_vec) * params[0]

    if mode == "phi":
        k0 = np.ones_like(x_vec) * np.tan(params[1])
    elif mode == "k":
        k0 = np.ones_like(x_vec) * params[1]
    else:
        assert False, f"Unknown mode '{mode}'"

    return y0 + k0 * x_vec


def straight_line_propagator_stepwise_2D_scatter_yphi(params, geo_pos, is_scatter):
    return straight_line_propagator_stepwise_2D_scatter(
        params, geo_pos, is_scatter, "phi"
    )


def straight_line_propagator_stepwise_2D_scatter_yk(params, geo_pos, is_scatter):
    return straight_line_propagator_stepwise_2D_scatter(
        params, geo_pos, is_scatter, "k"
    )


# This one works for y-phi and y-k descriptions
# The y-k-description will be converted to y-phi internally
def straight_line_propagator_stepwise_2D_scatter(
    params, geo_pos, is_scatter, mode="phi"
):
    assert len(geo_pos) == len(is_scatter)

    if mode == "phi":
        current_x_y_phi = [0, params[0], params[1]]
    elif mode == "k":
        current_x_y_phi = [0, params[0], np.arctan(params[1])]
    else:
        assert False, f"Unknown mode '{mode}'"

    new_x_y_phi = current_x_y_phi.copy()
    scatter_params = params[2:]

    y_vec = np.zeros_like(geo_pos)
    i_s = 0  # to count through the scattering surfaces
    for g in range(len(geo_pos)):
        # make hit
        new_x_y_phi = current_x_y_phi.copy()
        new_x_y_phi[0] = geo_pos[g]

        dx = new_x_y_phi[0] - current_x_y_phi[0]
        new_x_y_phi[1] = current_x_y_phi[1] + np.tan(current_x_y_phi[2]) * dx

        if is_scatter[g]:
            new_x_y_phi[2] = current_x_y_phi[2] + scatter_params[i_s]
            i_s += 1
        else:
            new_x_y_phi[2] = current_x_y_phi[2]

        y_vec[g] = new_x_y_phi[1]
        current_x_y_phi = new_x_y_phi.copy()

    return y_vec
