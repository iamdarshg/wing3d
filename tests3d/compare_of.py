"""OpenFOAM vs wing3d comparison report (forces + Cp cuts + verdicts)."""
import os
import sys
sys.path.insert(0, 'D:\\CodeProjects\\cfd')


def load_foam_log(logpath, last=100):
    from wing3d.openfoam import parse_log_forces
    return parse_log_forces(logpath, last=last)


def report_wing(foam, panel_inv, panel_visc):
    print('=== NACA0012 AR6, alpha=4deg, Re 3e6 (SST) ===')
    print(f'OpenFOAM RANS : CL={foam["CL"]:.4f} (+/-{foam["CL_std"]:.4f}) '
          f'CD={foam["CD"]:.5f} (p {foam["CDp"]:.5f} + v {foam["CDv"]:.5f})')
    print(f'wing3d invisc: CL={panel_inv["CL"]:.4f} CDp={panel_inv["CDp"]:.5f}')
    print(f'wing3d coupled: CL={panel_visc["CL"]:.4f} '
          f'CD={panel_visc["CD"]:.5f}')
    dcl = 100 * (panel_visc['CL'] - foam['CL']) / foam['CL']
    dcd = 100 * (panel_visc['CD'] - foam['CD']) / foam['CD']
    print(f'delta coupled vs RANS: CL {dcl:+.1f}%, CD {dcd:+.1f}%')
    for name, ok, detail in verdicts_wing(foam, panel_visc):
        print(('PASS' if ok else 'FAIL'), name, '-', detail)


def verdicts_wing(foam, panel):
    out = []
    dcl = abs(panel['CL'] - foam['CL']) / max(abs(foam['CL']), 1e-9)
    out.append(('CL within 15% of RANS', dcl < 0.15,
                f'dCL={dcl * 100:.1f}%'))
    dcd = abs(panel['CD'] - foam['CD']) / max(abs(foam['CD']), 1e-9)
    out.append(('CD within 30% of RANS', dcd < 0.30,
                f'dCD={dcd * 100:.1f}%'))
    out.append(('RANS lift positive & sane', foam['CL'] > 0.2,
                f"CL={foam['CL']:.3f}"))
    return out


def report_car(foam, panel):
    print('=== Ahmed-like slantback, alpha=0 (SST) ===')
    print(f'OpenFOAM RANS : CD={foam["CD"]:.4f} '
          f'(p {foam["CDp"]:.4f} + v {foam["CDv"]:.4f})')
    print(f'wing3d invisc: CDp={panel["CDp"]:.4f}')
    print('(inviscid omits separation: expect large gap on bluff bodies; '
          'RANS is the reference)')


if __name__ == '__main__':
    import json
    base = 'D:\\CodeProjects\\cfd\\cases\\openfoam'
    # filled by the runner once logs exist
    print('(run with parsed results; see compare runs)')
