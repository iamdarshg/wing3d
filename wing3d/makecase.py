"""Build a verified central-Euler case from proven templates + writer 0/ files.

Protocol (post 2026-10-08 audit): NEVER hand-sed BCs. Generate,
of_verify, deploy, verify-in-container, launch.
"""
import os
import shutil
import numpy as np

TEMPL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                     'cases', 'openfoam', 'templates', 'central')


def make_euler_central(root, mach, alpha_deg, p_inf=1e5, T_inf=300.0,
                       span_in=1.2, end_time=0.3):
    """Overwrite system/constant/0 of an existing case dir (mesh kept).

    Proven templates (biconic fvSchemes/fvSolution, LTS controlDict,
    janaf thermo) + freestream U at (mach, alpha). Slip wall, laminar.
    NOTE: template thermo is N2-janaf (R=296.8, from the tutorial) so
    U/rhoInf use R=296.8 consistently (~3% off air; documented, tiny
    vs mesh error).
    """
    Rgas = 8314.0 / 28.01348
    a_inf = float(np.sqrt(1.4 * Rgas * T_inf))
    Umag = float(mach * a_inf)
    a = float(np.radians(alpha_deg))
    ux, uz = Umag * np.cos(a), Umag * np.sin(a)
    rho_inf = p_inf / (Rgas * T_inf)
    for f in ['fvSchemes', 'fvSolution', 'controlDict']:
        shutil.copy(os.path.join(TEMPL, f),
                    os.path.join(root, 'system', f))
    shutil.copy(os.path.join(TEMPL, 'thermophysicalProperties'),
                os.path.join(root, 'constant', 'thermophysicalProperties'))
    with open(os.path.join(root, 'constant', 'turbulenceProperties'),
              'w', newline='\n') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; '
                'object turbulenceProperties; }\nsimulationType laminar;\n')
    with open(os.path.join(root, '0', 'U'), 'w', newline='\n') as f:
        f.write('FoamFile { version 2.0; format ascii; '
                'class volVectorField; object U; }\n'
                'dimensions [0 1 -1 0 0 0 0];\n'
                'internalField uniform (%.3f 0 %.3f);\n'
                'boundaryField {\n'
                ' farfield { type freestream; freestreamValue uniform '
                '(%.3f 0 %.3f); }\n'
                ' symm { type symmetry; }\n'
                ' "slab.stl" { type slip; }\n}\n' % (ux, uz, ux, uz))
    with open(os.path.join(root, '0', 'p'), 'w', newline='\n') as f:
        f.write('FoamFile { version 2.0; format ascii; '
                'class volScalarField; object p; }\n'
                'dimensions [1 -1 -2 0 0 0 0];\n'
                'internalField uniform %.0f;\n'
                'boundaryField {\n'
                ' farfield { type freestreamPressure;\n'
                '  freestreamValue uniform %.0f; }\n'
                ' symm { type symmetry; }\n'
                ' "slab.stl" { type zeroGradient; }\n}\n'
                % (p_inf, p_inf))
    with open(os.path.join(root, '0', 'T'), 'w', newline='\n') as f:
        f.write('FoamFile { version 2.0; format ascii; '
                'class volScalarField; object T; }\n'
                'dimensions [0 0 0 1 0 0 0];\n'
                'internalField uniform %.0f;\n'
                'boundaryField {\n'
                ' farfield { type freestream; freestreamValue uniform '
                '%.0f; }\n'
                ' symm { type symmetry; }\n'
                ' "slab.stl" { type zeroGradient; }\n}\n' % (T_inf, T_inf))
    for gone in ['nut', 'nuTilda', 'alphat']:
        p = os.path.join(root, '0', gone)
        if os.path.exists(p):
            os.remove(p)
    # controlDict: LTS endTime + forces numbers
    cd = os.path.join(root, 'system', 'controlDict')
    txt = open(cd).read()
    import re
    txt = re.sub(r'stopAt endTime; endTime [^;]+;',
                 'stopAt endTime; endTime %s;' % end_time, txt)
    txt = re.sub(r'magUInf [\d.eE+-]+', 'magUInf %.3f' % Umag, txt)
    txt = re.sub(r'rhoInf [\d.eE+-]+', 'rhoInf %.4f' % rho_inf, txt)
    txt = re.sub(r'(deltaT\s+)[\d.eE+-]+', r'\g<1>2e-5', txt)
    txt = re.sub(r'maxCo\s+[\d.eE+-]+', 'maxCo 0.5', txt)
    txt = re.sub(r'maxDeltaT\s+[\d.eE+-]+', 'maxDeltaT 1e-4', txt)
    txt = re.sub(r'startFrom\s+\w+;', 'startFrom startTime;', txt)
    open(cd, 'w', newline='\n').write(txt)
    print('central case %s: M=%.2f a=%.2f U=%.1f' % (root, mach,
                                                     alpha_deg, Umag))
    return {'U': Umag, 'rho': rho_inf}


if __name__ == '__main__':
    import sys
    make_euler_central(sys.argv[1], float(sys.argv[2]),
                       float(sys.argv[3]))
