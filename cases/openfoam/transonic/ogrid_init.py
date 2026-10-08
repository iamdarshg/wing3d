import numpy as np

Rgas = 8314.0 / 28.01348
T_inf = 300.0
p_inf = 1e5
a_inf = float(np.sqrt(1.4 * Rgas * T_inf))
MACH, ALPHA = 0.8, 1.25
Umag = MACH * a_inf
a = np.radians(ALPHA)
ux, uz = Umag * np.cos(a), Umag * np.sin(a)
D = 'D:/CodeProjects/cfd/cases/openfoam/ogrid0'
import os
os.makedirs(D, exist_ok=True)
open(D + '/U', 'w', newline='\n').write(
    'FoamFile { version 2.0; format ascii; class volVectorField; '
    'object U; }\ndimensions [0 1 -1 0 0 0 0];\n'
    'internalField uniform (%.3f 0 %.3f);\nboundaryField {\n'
    ' farfield { type freestream; freestreamValue uniform (%.3f 0 %.3f); }\n'
    ' airfoil { type slip; }\n frontback { type empty; }\n}\n'
    % (ux, uz, ux, uz))
open(D + '/p', 'w', newline='\n').write(
    'FoamFile { version 2.0; format ascii; class volScalarField; '
    'object p; }\ndimensions [1 -1 -2 0 0 0 0];\n'
    'internalField uniform %.0f;\nboundaryField {\n'
    ' farfield { type freestreamPressure; freestreamValue uniform %.0f; }\n'
    ' airfoil { type zeroGradient; }\n frontback { type empty; }\n}\n'
    % (p_inf, p_inf))
open(D + '/T', 'w', newline='\n').write(
    'FoamFile { version 2.0; format ascii; class volScalarField; '
    'object T; }\ndimensions [0 0 0 1 0 0 0];\n'
    'internalField uniform %.0f;\nboundaryField {\n'
    ' farfield { type freestream; freestreamValue uniform %.0f; }\n'
    ' airfoil { type zeroGradient; }\n frontback { type empty; }\n}\n'
    % (T_inf, T_inf))
print('U=%.1f rho=%.4f' % (Umag, p_inf / (Rgas * T_inf)))
