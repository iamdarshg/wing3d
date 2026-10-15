"""Quick TSD baseline probe: M, span, free_k flag -> CL at campaign grid."""
import sys
import time
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_m_sup import run_sup
from tsd_campaign import build_profile

M = float(sys.argv[1])
span = float(sys.argv[2]) if len(sys.argv) > 2 else 1.6
fk = len(sys.argv) > 3 and sys.argv[3] == 'fk'
prof = build_profile('naca0012')
t0 = time.time()
L, xc, ix, j0 = run_sup(M, span=span, epochs=4, free_k=fk, profile=prof)
print('PROBE M=%.3f span=%.1f free_k=%s CL=%.4f el=%.0fs' % (
    M, span, fk, L['CL'], time.time() - t0))
