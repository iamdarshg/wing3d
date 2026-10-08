"""C kernel access (cffi ABI) with numpy fallback.

Build once: python -m wing3d.csrc.build  (needs gcc/cl, ships one DLL).
Disable: set WING3D_NO_C=1.

NOTE: converted arrays are bound to locals before .ctypes.data --
inline temporaries would be freed before the C call runs (dangling
pointers). Do not "simplify" back to inline conversion.
"""
import os

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc')
LIB = os.path.join(HERE, 'panel_rows' + ('.dll' if os.name == 'nt'
                                         else '.so'))
_C_OK = False
_lib = None

if not os.environ.get('WING3D_NO_C'):
    try:
        from cffi import FFI as _FFI
        _ffi = _FFI()
        _ffi.cdef("""
void solid_angle_block(const double *P, int np_, const double *fan,
    const unsigned char *fmask, int n, int F, double *om);
void source_pot_block(const double *P, int np_, const double *qp,
    const double *qw, int n, int Q, const double *area,
    const double *size, const double *cent, double *out);
void source_vel_block(const double *P, int np_, const double *qp,
    const double *qw, int n, int Q, const double *area,
    const double *size, const double *cent, double *out);
void loop_vel_block(const double *P, int np_, const double *ed,
    const unsigned char *emask, int n, int E, double *out);
""")
        if os.path.exists(LIB):
            _lib = _ffi.dlopen(LIB)
            _C_OK = True
    except Exception:
        _C_OK = False


def _req(a, dtype=None):
    import numpy as np
    return np.ascontiguousarray(a, dtype=dtype or np.float64)


def solid_angle_block(P, fan, fmask):
    """om[t, j] solid angle (unscaled sum). Numpy fallback if no C."""
    import numpy as np
    Pc = _req(P)
    fanc = _req(fan)
    fmc = np.ascontiguousarray(fmask, dtype=np.uint8)
    n = fanc.shape[0]
    F = fanc.shape[1]
    if not _C_OK:
        from .solver import _solid_angle_rows
        return np.array([_solid_angle_rows(p, fan, fmask) for p in Pc])
    om = np.zeros((len(Pc), n))
    _lib.solid_angle_block(
        _ffi.cast('const double *', Pc.ctypes.data),
        len(Pc),
        _ffi.cast('const double *', fanc.ctypes.data),
        _ffi.cast('const unsigned char *', fmc.ctypes.data),
        n, F,
        _ffi.cast('double *', om.ctypes.data))
    return om


def source_pot_block(P, qp, qw, area, size, cent):
    import numpy as np
    Pc = _req(P)
    qpc = _req(qp)
    qwc = _req(qw)
    areac = _req(area)
    sizec = _req(size)
    centc = _req(cent)
    n = qpc.shape[0]
    if not _C_OK:
        from .solver import _source_pot_rows
        return np.array([_source_pot_rows(p, qp, qw, area, size)
                         for p in Pc])
    out = np.zeros((len(Pc), n))
    _lib.source_pot_block(
        _ffi.cast('const double *', Pc.ctypes.data), len(Pc),
        _ffi.cast('const double *', qpc.ctypes.data),
        _ffi.cast('const double *', qwc.ctypes.data), n, qpc.shape[1],
        _ffi.cast('const double *', areac.ctypes.data),
        _ffi.cast('const double *', sizec.ctypes.data),
        _ffi.cast('const double *', centc.ctypes.data),
        _ffi.cast('double *', out.ctypes.data))
    return out


def source_vel_block(P, qp, qw, area, size, cent):
    import numpy as np
    Pc = _req(P)
    qpc = _req(qp)
    qwc = _req(qw)
    areac = _req(area)
    sizec = _req(size)
    centc = _req(cent)
    n = qpc.shape[0]
    if not _C_OK:
        from .solver import _source_vel_rows
        return np.array([_source_vel_rows(p, qp, qw, area, size)
                         for p in Pc])
    out = np.zeros((len(Pc), n, 3))
    _lib.source_vel_block(
        _ffi.cast('const double *', Pc.ctypes.data), len(Pc),
        _ffi.cast('const double *', qpc.ctypes.data),
        _ffi.cast('const double *', qwc.ctypes.data), n, qpc.shape[1],
        _ffi.cast('const double *', areac.ctypes.data),
        _ffi.cast('const double *', sizec.ctypes.data),
        _ffi.cast('const double *', centc.ctypes.data),
        _ffi.cast('double *', out.ctypes.data))
    return out


def loop_vel_block(P, ed, emask):
    import numpy as np
    Pc = _req(P)
    edc = _req(ed)
    emc = np.ascontiguousarray(emask, dtype=np.uint8)
    n = edc.shape[0]
    if not _C_OK:
        from .solver import _loop_vel_rows
        return np.array([_loop_vel_rows(p, ed, emask) for p in Pc])
    out = np.zeros((len(Pc), n, 3))
    _lib.loop_vel_block(
        _ffi.cast('const double *', Pc.ctypes.data), len(Pc),
        _ffi.cast('const double *', edc.ctypes.data),
        _ffi.cast('const unsigned char *', emc.ctypes.data),
        n, edc.shape[1],
        _ffi.cast('double *', out.ctypes.data))
    return out
