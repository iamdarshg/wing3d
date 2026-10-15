"""Compare CLI-extracted RAE section vs the campaign table section."""
import sys
import numpy as np
sys.path.insert(0, r'D:\CodeProjects\cfd')
sys.path.insert(0, r'D:\CodeProjects\cfd\tests3d')

from wing3d.tsdcli import extract_section
from wing3d.geometry import read_stl, section_loop
from gen_transonic_stl import rae2822_xy


def main():
    mesh = read_stl(r'D:\CodeProjects\cfd\cases\openfoam\transonic\stl\rae2822.stl')
    sec = extract_section(mesh)

    xu, yu, xl, yl = rae2822_xy()
    if xu[0] > xu[-1]:
        xu, yu = xu[::-1], yu[::-1]
    if xl[0] > xl[-1]:
        xl, yl = xl[::-1], yl[::-1]

    print('table x range [%.4f, %.4f]  y range [%.4f, %.4f]'
          % (xu[0], xu[-1], min(yu.min(), yl.min()), max(yu.max(), yl.max())))
    print('cli   x range [%.4f, %.4f]  y range [%.4f, %.4f]'
          % (sec['x'][0], sec['x'][-1], min(sec['yu'].min(), sec['yl'].min()),
             max(sec['yu'].max(), sec['yl'].max())))

    st = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.76, 0.8, 0.9,
          0.95, 0.98, 0.995]
    print('%-6s %9s %9s %8s | %9s %9s %8s | %8s %8s' %
          ('x', 'tbl_u', 'cli_u', 'dU', 'tbl_l', 'cli_l', 'dL',
           'tbl_cb', 'cli_cb'))
    for s in st:
        tu = np.interp(s, xu, yu)
        tl = np.interp(s, xl, yl)
        cu = np.interp(s, sec['x'], sec['yu'])
        cl = np.interp(s, sec['x'], sec['yl'])
        print('%-6.3f %+9.4f %+9.4f %+8.4f | %+9.4f %+9.4f %+8.4f | %+8.4f %+8.4f'
              % (s, tu, cu, cu - tu, tl, cl, cl - tl,
                 0.5 * (tu + tl), 0.5 * (cu + cl)))

    loop_a = section_loop(*rae2822_xy(), 24)
    loop_b = section_loop(sec['x'], sec['yu'], sec['x'], sec['yl'], 24)
    print('campaign loop pts', loop_a.shape, ' cli loop pts', loop_b.shape)
    d = np.hypot(loop_a[:, 0] - loop_b[:, 0], loop_a[:, 1] - loop_b[:, 1])
    print('loop max |diff|=%.5f  mean |diff|=%.5f' % (d.max(), d.mean()))
    for i in range(0, len(loop_a), 3):
        print('  pt %2d  camp(%.4f,%+.4f)  cli(%.4f,%+.4f)  d=%.4f'
              % (i, loop_a[i, 0], loop_a[i, 1], loop_b[i, 0], loop_b[i, 1],
                 d[i]))


if __name__ == '__main__':
    main()
