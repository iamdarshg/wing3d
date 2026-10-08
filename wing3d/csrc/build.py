"""Build panel_rows DLL with gcc (Windows) or cc (posix)."""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'panel_rows.c')
OUT = os.path.join(HERE, 'panel_rows' + ('.dll' if sys.platform == 'win32'
                                         else '.so'))


def build():
    cc = shutil.which('gcc') or shutil.which('cc') or shutil.which('cl')
    if cc is None:
        raise RuntimeError('no C compiler found (gcc/cc/cl)')
    base = os.path.basename(cc).lower()
    if 'cl' in base:
        cmd = [cc, '/O2', '/LD', SRC, '/Fe' + OUT]
    else:
        cmd = [cc, '-O2', '-shared', '-o', OUT, SRC, '-lm']
    print('building:', ' '.join(cmd))
    subprocess.check_call(cmd)
    print('wrote', OUT)
    return OUT


if __name__ == '__main__':
    build()
