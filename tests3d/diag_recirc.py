"""Run the recirc (TSFOIL-style) path with timing; CLI: M [nl] [epochs]."""
import sys
import time
sys.path.insert(0, 'D:\\CodeProjects\\cfd')
sys.path.insert(0, 'D:\\CodeProjects\\cfd\\tests3d')
from tsd_recirc import run_recirc

M = float(sys.argv[1]) if len(sys.argv) > 1 else 0.85
nl = len(sys.argv) > 2 and sys.argv[2] == 'nl'
epochs = int(sys.argv[3]) if len(sys.argv) > 3 else 6
t0 = time.time()
s = run_recirc(M, alpha_deg=1.25, linear=not nl, epochs=epochs)
import numpy as np
G = s.wake_G[s.wing_jy]
print('final: max|G|=%.4f  elapsed=%.1f s' % (float(np.abs(G).max()),
                                              time.time() - t0))
