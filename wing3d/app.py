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


def build_case(kind, alpha=4.0, n_chord=20, n_span=10, stl_bytes=None):
    import numpy as np
    from wing3d.geometry import build_wing, read_stl
    from wing3d.solver import build_wake_from_meta
    from wing3d.primitives import box
    from wing3d.shapes import build_car
    if kind == 'wing':
        mesh = build_wing('0012', span=6.0, chord=1.0, n_chord=n_chord,
                          n_span=n_span)
        wakes = build_wake_from_meta(mesh)
        sref, strips = 6.0, True
    elif kind == 'car':
        mesh = build_car()
        wakes, sref, strips = [], 0.389 * 0.288, False
    elif kind == 'cube':
        mesh = box(size=(1, 1, 1), n=(10, 10, 10))
        wakes, sref, strips = [], 1.0, False
    elif kind == 'stl':
        import tempfile
        fd, path = tempfile.mkstemp(suffix='.stl')
        os.write(fd, stl_bytes)
        os.close(fd)
        mesh = read_stl(path)
        os.remove(path)
        from wing3d.geometry import decimate
        if mesh.npanels > 2500:
            mesh = decimate(mesh, target=2000)
        wakes, sref, strips = [], 1.0, False
    else:
        raise ValueError('unknown geometry')
    a = np.radians(alpha)
    vinf = [np.cos(a), 0, np.sin(a)]
    return mesh, wakes, vinf, sref, strips


def solve_case(mesh, wakes, vinf, sref):
    import numpy as np
    from wing3d.solver import solve, assemble
    from wing3d.forces import pressure_forces
    A, Brow = assemble(mesh, wakes,
                       kutta_mode='doublet' if wakes else 'doublet')
    res = solve(mesh, wakes, vinf, A=A, Brow=Brow, kutta_mode='doublet')
    f = pressure_forces(mesh, res['cp'], vinf, sref=sref)
    return res, f


PAGE = """<html><head><title>wing3d</title></head><body>
<h1>wing3d — full-3D panel + viscous simulator</h1>
<form action="/run" method="post" enctype="multipart/form-data">
Geometry: <select name="kind">
<option value="wing">NACA0012 wing AR6</option>
<option value="car">Ahmed-like car</option>
<option value="cube">cube</option>
<option value="stl">STL upload</option>
</select><br>
Alpha (deg): <input name="alpha" value="4" size="5"><br>
Chord panels: <input name="n_chord" value="20" size="5">
Span panels: <input name="n_span" value="10" size="5"><br>
STL file: <input type="file" name="stl"><br>
<input type="submit" value="Run">
</form>
<p>Validation: sphere/cube/ellipsoid analytic + wing polar in
tests3d/test_shapes.py. OpenFOAM comparison: cases/openfoam/.</p>
</body></html>"""

RESULT = """<html><head><title>wing3d result</title></head><body>
<h1>Result: {kind} at alpha={alpha}</h1>
<p>CL={CL:.4f} CDp={CDp:.5f} CS={CS:.4f} (sref={sref}) panels={npan}</p>
<img src="/img/{a}"><br>
<img src="/img/{b}"><br>
<a href="/vtk/{v}">Download VTK (ParaView)</a><br>
<a href="/">Back</a></body></html>"""


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/':
            self._html(PAGE)
        elif self.path.startswith('/img/'):
            name = self.path[5:]
            self._file(os.path.join(OUT, name), 'image/png')
        elif self.path.startswith('/vtk/'):
            name = self.path[5:]
            self._file(os.path.join(OUT, name), 'application/octet-stream')
        else:
            self.send_error(404)

    def do_POST(self):
        import cgi
        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers,
                                environ={'REQUEST_METHOD': 'POST'})
        kind = form.getvalue('kind', 'wing')
        alpha = float(form.getvalue('alpha', '4'))
        nc = int(form.getvalue('n_chord', '20'))
        ns = int(form.getvalue('n_span', '10'))
        stl_bytes = None
        if 'stl' in form and getattr(form['stl'], 'file', None):
            stl_bytes = form['stl'].file.read()
            kind = 'stl'
        mesh, wakes, vinf, sref, _ = build_case(kind, alpha, nc, ns,
                                                stl_bytes)
        res, f = solve_case(mesh, wakes, vinf, sref)
        from wing3d.viz import plot_surface_cp, write_vtk
        tag = f'{kind}_a{alpha}'.replace('.', 'p')
        a = tag + '_cp.png'
        plot_surface_cp(mesh, res['cp'], os.path.join(OUT, a),
                        title=f'{kind} Cp alpha={alpha}')
        v = tag + '.vtk'
        write_vtk(mesh, res['cp'], os.path.join(OUT, v))
        self._html(RESULT.format(kind=kind, alpha=alpha, CL=f['CL'],
                                 CDp=f['CDp'], CS=f['CS'], sref=sref,
                                 npan=mesh.npanels, a=a, b=a, v=v))

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
