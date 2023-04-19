# Alexander J. Pfleger
# 2023-03-21
#
# SymPy sheet to check all derivations

from sympy import *

x, xs, y, phi, theta, y0, phi0, theta0, y1, phi1 = symbols(
    "x xs y phi theta y0 phi0 theta0 y1 phi1"
)

init_printing(use_unicode=True)

# before scattering
y = y0 + tan(phi0) * x
phi = phi0
theta = theta0

diff(y, y0)
simplify(diff(y, phi0))
diff(y, theta0)


# after scattering
phi1 = phi0 + theta0
# y1 + tan(phi1) * xs =  y0 + tan(phi0) * xs
y1 = -tan(phi1) * xs + y0 + tan(phi0) * xs

y = y1 + tan(phi1) * x
phi = phi1
theta = theta0

diff(y, y0)
diff(y, phi0)
simplify(diff(y, theta0))


simplify(diff(y, phi0) * diff(y, theta0))
