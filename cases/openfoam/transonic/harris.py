"""Harris NACA0012 data (NASA TM-81927, via published tables).

M=0.749, alpha=1.99deg, Re=2.0e6 (McDevitt dataset, OCR'd from the
report text -- values cross-checked for monotonicity; +-0.01 caveat).
Upper: shock between x=0.50 (Cp -1.187) and x=0.55 (Cp -0.340).
Source: NASA-TM-81927 / AGARD propagation; use as transonic truth
for pocket/shock/CL (integrate below).
"""
import numpy as np

# (x/c, Cp) upper surface, M0.749 a1.99
HARRIS_UPPER = np.array([
    [0.000, 1.081], [0.025, -0.507], [0.050, -0.721], [0.075, -0.881],
    [0.100, -0.956], [0.150, -1.064], [0.200, -1.135], [0.250, -1.175],
    [0.300, -1.211], [0.350, -1.230], [0.400, -1.220], [0.450, -1.210],
    [0.500, -1.187], [0.550, -0.340], [0.600, -0.219], [0.650, -0.165],
    [0.700, -0.126], [0.750, -0.093], [0.800, -0.049], [0.850, -0.006],
    [0.900, 0.049], [0.925, 0.082], [0.950, 0.115], [0.975, 0.169],
    [1.000, 0.222],
])
# (x/c, Cp) lower surface (partial stations from same table)
HARRIS_LOWER = np.array([
    [0.000, 1.081], [0.025, 0.172], [0.050, -0.087], [0.100, -0.335],
    [0.200, -0.476], [0.300, -0.468], [0.400, -0.380], [0.500, -0.305],
    [0.600, -0.220], [0.700, -0.149], [0.800, -0.072], [0.850, -0.020],
    [0.900, 0.036], [0.950, 0.101], [1.000, 0.174],
])


def integrate_cl_cd(xu, cpu, xl, cpl):
    """Section cl (trapezoid of dCp) + crude pressure-drag estimate."""
    # order by x
    ou = np.argsort(xu)
    ol = np.argsort(xl)
    xu, cpu = xu[ou], cpu[ou]
    xl, cpl = xl[ol], cpl[ol]
    # common grid
    xg = np.linspace(0, 1, 401)
    cpu_i = np.interp(xg, xu, cpu)
    cpl_i = np.interp(xg, xl, cpl)
    cl = float(np.trapezoid(cpl_i - cpu_i, xg))
    return cl


if __name__ == '__main__':
    cl = integrate_cl_cd(HARRIS_UPPER[:, 0], HARRIS_UPPER[:, 1],
                         HARRIS_LOWER[:, 0], HARRIS_LOWER[:, 1])
    print('Harris M0.749 a1.99: cl = %.4f' % cl)
    print('shock between x=0.50 and 0.55; Cpmin = %.2f' % HARRIS_UPPER[:, 1].min())
