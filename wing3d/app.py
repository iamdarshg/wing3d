"""wing3d local web UI (stdlib only, offline). Run: python -m wing3d.app."""
import io
import os
import sys
import urllib.parse
from http.server import BaseHTTPRequestHandler, HTTPServer

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'cases', 'ui_out')
os.makedirs(OUT, exist_ok=True)

CSS = """
:root{color-scheme:dark;--bg:#0c1117;--surface:#111820;--panel:#10171e;
--border:#26303a;--text:#e4ecef;--muted:#8796a5;--accent:#89e5d0;
--mono:"SFMono-Regular",Consolas,"Liberation Mono",monospace;
font-family:Inter,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
font-size:13px}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text)}
a{color:var(--accent);text-decoration:none}
.topbar{height:72px;border-bottom:1px solid var(--border);display:flex;
align-items:center;justify-content:space-between;padding:0 30px;background:#0e151c}
.brand{display:flex;gap:11px;align-items:center;color:var(--text);
font-size:23px;font-weight:600;letter-spacing:-1px}
.brand small{display:block;font-size:8px;font-weight:500;letter-spacing:1.5px;
margin-top:4px;color:#97a7b3}
.tag{font-family:var(--mono);font-size:9px;letter-spacing:.5px;padding:4px 6px;
border:1px solid #2a3c41;border-radius:4px;color:#9aafb4;white-space:nowrap}
.app-shell{display:grid;grid-template-columns:282px minmax(0,1fr)}
.sidebar{padding:26px 22px 20px;border-right:1px solid var(--border);background:#0e151c}
.sidebar h1{font-size:19px;letter-spacing:-.5px;font-weight:550;margin:15px 0 6px}
.field-label{font-size:10px;color:#99aab6;display:block;margin-bottom:7px}
.section-label{display:flex;align-items:center;justify-content:space-between;
margin:25px 0 12px;font-size:10px;letter-spacing:1.2px;color:#9baab7}
input,select{min-width:0;border:1px solid var(--border);border-radius:6px;
background:#0b1219;color:var(--text);outline:none;padding:9px 10px;width:100%}
input:focus,select:focus{border-color:var(--accent)}
input[type=number]{font-family:var(--mono);font-size:12px}
.quiet{border:1px solid var(--border);border-radius:6px;padding:8px 12px;
font-size:11px;background:#121b23;color:var(--text);cursor:pointer}
.quiet:hover{border-color:#627e89}
main{padding:28px 30px 16px;min-width:0;max-width:1800px;width:100%;margin:auto}
.workspace-heading h2{font-size:25px;font-weight:500;letter-spacing:-.65px;margin:10px 0 0}
.status-line{display:flex;align-items:center;gap:11px;margin:18px 0 20px;
font-size:10px;color:var(--muted);min-height:22px}
.status-pill{font-family:var(--mono);color:#a6d6c8;background:#17302c;
border:1px solid #2c4a41;border-radius:4px;padding:4px 7px;font-size:9px;
text-transform:uppercase;letter-spacing:.6px}
.metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-bottom:22px}
.metric{background:linear-gradient(130deg,#151e27,#111820);border:1px solid var(--border);
border-radius:8px;padding:15px 17px}
.metric>span{font-size:10px;color:#9daeba;display:flex;justify-content:space-between}
.metric strong{display:block;font-family:var(--mono);font-size:29px;font-weight:400;
letter-spacing:-1px;color:var(--accent);margin:12px 0 7px}
.metric small{font-size:9px;color:#778c9b}
.plot-panel{border:1px solid var(--border);border-radius:8px;background:var(--panel);
overflow:hidden;margin-bottom:18px}
.panel-heading{height:51px;display:flex;align-items:center;
justify-content:space-between;padding:0 18px;border-bottom:1px solid #26303a88}
.panel-heading h3{font-size:12px;font-weight:500;margin:0}
.panel-index{font-family:var(--mono);font-size:10px;color:#678591}
.plot-panel img{width:100%;display:block}
.subtle{color:var(--muted);font-size:11px}
"""

PAGE = """<html><head><title>wing3d · 3D workbench</title><meta name="theme-color" content="#0c1117">
<style>{css}</style></head><body>
<header class="topbar"><div class="brand">wing3d<small>FULL-3D PANEL + VISCOUS SIMULATOR</small></div>
<div><span class="tag">LOCAL</span> <span class="tag">NO UPLOAD</span></div></header>
<div class="app-shell">
<aside class="sidebar"><h1>Test case</h1>
<form action="/run" method="post" enctype="multipart/form-data">
<div class="section-label"><span>GEOMETRY</span></div>
<label class="field-label" for="kind">Configuration</label>
<select name="kind" id="kind">
<option value="wing">NACA0012 wing AR6</option>
<option value="car">Ahmed-like car</option>
<option value="cube">cube</option>
<option value="sphere">sphere</option>
<option value="stl">STL upload</option>
</select>
<div class="section-label"><span>FLOW CONDITIONS</span><span class="tag">INCOMPRESSIBLE</span></div>
<label class="field-label" for="alpha">Angle of attack (deg)</label>
<input name="alpha" id="alpha" type="number" value="4" step="any">
<label class="field-label" for="rey">Reynolds number</label>
<input name="rey" id="rey" type="number" value="3000000" step="any">
<div class="section-label"><span>MESH</span></div>
<label class="field-label" for="n_chord">Chord panels</label>
<input name="n_chord" id="n_chord" type="number" value="20">
<label class="field-label" for="n_span">Span panels</label>
<input name="n_span" id="n_span" type="number" value="10">
<div class="section-label"><span>STL UPLOAD</span></div>
<input type="file" name="stl">
<p><button class="quiet" type="submit">Run analysis</button></p>
</form>
<p class="subtle">Validation: sphere/cube/ellipsoid analytic + wing polar in
tests3d/test_shapes.py. OpenFOAM comparison: cases/openfoam/.</p>
</aside>
<main><div class="workspace-heading"><h2>3D visualization</h2></div>
<div class="status-line"><span class="status-pill">ready</span><span>pick a case on the left</span></div>
<div class="metrics">
<div class="metric"><span>CL</span><strong>—</strong><small>lift coefficient</small></div>
<div class="metric"><span>CD</span><strong>—</strong><small>drag coefficient</small></div>
<div class="metric"><span>L/D</span><strong>—</strong><small>aerodynamic efficiency</small></div>
<div class="metric"><span>PANELS</span><strong>—</strong><small>surface mesh</small></div>
</div></main></div></body></html>"""

RESULT = """<html><head><title>wing3d result</title><meta name="theme-color" content="#0c1117">
<style>{css}</style></head><body>
<header class="topbar"><div class="brand">wing3d<small>FULL-3D PANEL + VISCOUS SIMULATOR</small></div>
<div><a href="/">← New case</a></div></header>
<div class="app-shell">
<aside class="sidebar"><h1>Case summary</h1>
<div class="section-label"><span>CONFIGURATION</span></div>
<p class="subtle">Geometry: {kind}<br>Alpha: {alpha}°<br>Re: {rey}<br>Panels: {npan}</p>
<div class="section-label"><span>EXPORT</span></div>
<p><a class="quiet" href="/vtk/{v}">↓ Download VTK (ParaView)</a></p>
</aside>
<main><div class="workspace-heading"><h2>3D visualization</h2></div>
<div class="status-line"><span class="status-pill">solved</span><span>{elapsed}</span></div>
<div class="metrics">
<div class="metric"><span>CL</span><strong>{CL:.4f}</strong><small>lift coefficient</small></div>
<div class="metric"><span>CD</span><strong>{CD:.5f}</strong><small>drag coefficient</small></div>
<div class="metric"><span>L/D</span><strong>{LD:.1f}</strong><small>aerodynamic efficiency</small></div>
<div class="metric"><span>PANELS</span><strong>{npan}</strong><small>surface mesh</small></div>
</div>
<div class="plot-panel"><div class="panel-heading"><h3>Surface pressure</h3>
<span class="panel-index">01 / Cp</span></div><img src="/img/{a}"></div>
<div class="plot-panel"><div class="panel-heading"><h3>Spanwise loading</h3>
<span class="panel-index">02 / cl(y)</span></div><img src="/img/{b}"></div>
</main></div></body></html>"""


def build_case(kind, alpha=4.0, n_chord=20, n_span=10, stl_bytes=None,
               rey=3e6):
    import numpy as np
    from wing3d.geometry import build_wing, read_stl
    from wing3d.solver import build_wake_from_meta
    from wing3d.primitives import box
    from wing3d.shapes import build_car
    from tests3d.analytic_sphere import sphere_mesh
    if kind == 'wing':
        mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=n_chord,
                          n_span=n_span)
        wakes = build_wake_from_meta(mesh)
        sref = 6.0
    elif kind == 'car':
        mesh = build_car()
        wakes, sref = [], 0.389 * 0.288
    elif kind == 'cube':
        mesh = box(size=(1, 1, 1), n=(10, 10, 10))
        wakes, sref = [], 1.0
    elif kind == 'sphere':
        mesh = sphere_mesh(12, 24)
        wakes, sref = [], np.pi
    elif kind == 'stl':
        import tempfile
        import os as _os
        fd, path = tempfile.mkstemp(suffix='.stl')
        _os.write(fd, stl_bytes)
        _os.close(fd)
        mesh = read_stl(path)
        _os.remove(path)
        from wing3d.geometry import decimate
        if mesh.npanels > 2500:
            mesh = decimate(mesh, target=2000)
        wakes, sref = [], 1.0
    else:
        raise ValueError('unknown geometry')
    a = np.radians(alpha)
    vinf = [np.cos(a), 0, np.sin(a)]
    return mesh, wakes, vinf, sref


def solve_case(mesh, wakes, vinf, sref):
    import time
    from wing3d.solver import solve, assemble
    from wing3d.forces import pressure_forces
    t0 = time.time()
    A, Brow = assemble(mesh, wakes,
                       kutta_mode='doublet' if wakes else 'doublet')
    res = solve(mesh, wakes, vinf, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], vinf, sref=sref)
    return res, f, time.time() - t0


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/':
            self._html(PAGE.format(css=CSS))
        elif self.path.startswith('/img/'):
            self._file(os.path.join(OUT, self.path[5:]), 'image/png')
        elif self.path.startswith('/vtk/'):
            self._file(os.path.join(OUT, self.path[5:]),
                       'application/octet-stream')
        else:
            self.send_error(404)

    def do_POST(self):
        import cgi
        import time
        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers,
                                environ={'REQUEST_METHOD': 'POST'})
        kind = form.getvalue('kind', 'wing')
        alpha = float(form.getvalue('alpha', '4'))
        rey = float(form.getvalue('rey', '3e6'))
        nc = int(form.getvalue('n_chord', '20'))
        ns = int(form.getvalue('n_span', '10'))
        stl_bytes = None
        if 'stl' in form and getattr(form['stl'], 'file', None):
            stl_bytes = form['stl'].file.read()
            kind = 'stl'
        t0 = time.time()
        mesh, wakes, vinf, sref = build_case(kind, alpha, nc, ns, stl_bytes,
                                             rey)
        res, f, _ = solve_case(mesh, wakes, vinf, sref)
        from wing3d.viz import plot_surface_cp, plot_spanload, write_vtk
        tag = f'{kind}_a{alpha}'.replace('.', 'p')
        a = tag + '_cp.png'
        plot_surface_cp(mesh, res['cp'], os.path.join(OUT, a),
                        title=f'{kind} Cp alpha={alpha}')
        import numpy as np
        C = mesh.centroid
        order = np.argsort(C[:, 1])
        b = tag + '_span.png'
        plot_spanload(C[order][:, 1], res['cp'][order],
                      os.path.join(OUT, b), title='Cp vs span')
        v = tag + '.vtk'
        write_vtk(mesh, res['cp'], os.path.join(OUT, v))
        ld = abs(f['CL'] / max(abs(f['CDp']), 1e-9))
        self._html(RESULT.format(css=CSS, kind=kind, alpha=alpha, rey=rey,
                                 CL=f['CL'], CD=f['CDp'], LD=ld,
                                 npan=mesh.npanels, a=a, b=b, v=v,
                                 elapsed=f'solved in {time.time()-t0:.1f}s'))

    def _html(self, s):
        b = s.encode()
        self.send_response(200)
        self.send_header('Content-Type', 'text/html')
        self.send_header('Content-Length', str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def _file(self, path, ctype):
        with open(path, 'rb') as fh:
            b = fh.read()
        self.send_response(200)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(b)))
        self.end_headers()
        self.wfile.write(b)


def main(port=8080):
    print(f'wing3d UI on http://localhost:{port}')
    HTTPServer(('127.0.0.1', port), Handler).serve_forever()


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8080)
