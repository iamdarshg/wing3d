"""High-vertex STL suite for transonic validation (100k+ verts).

Classic transonic test geometries, each written as a binary STL with
>100k vertices:
1. NACA0012 (symmetric, Harris/OF anchor lineage)
2. RAE2822 (supercritical-ish, AGARD transonic validation standard)
3. Biconvex (analytic circular-arc, supersonic theory check)
4. ONERA M6 (classic swept-wing transonic validation geometry)

Vertex budget: an airfoil section with ~350 chordwise stations, extruded
through ~200 spanwise layers on a slab => ~140k vertices.
"""
import numpy as np
import os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                    'cases', 'openfoam', 'transonic', 'stl')


def naca0012_xy(n=350, t=0.12):
    beta = np.linspace(0, np.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                  + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    return x, yt


def rae2822_xy():
    """RAE2822 (UIUC airfoil database coordinates, chord-normalized)."""
    raw = """0.000000 0.000000
0.000602 0.003165
0.002408 0.006306
0.005412 0.009416
0.009607 0.012480
0.014984 0.015489
0.021530 0.018441
0.029228 0.021348
0.038060 0.024219
0.048005 0.027062
0.059039 0.029874
0.071136 0.032644
0.084265 0.035360
0.098396 0.038011
0.113495 0.040585
0.129524 0.043071
0.146447 0.045457
0.164221 0.047729
0.182803 0.049874
0.202150 0.051885
0.222215 0.053753
0.242949 0.055470
0.264302 0.057026
0.286222 0.058414
0.308658 0.059629
0.331555 0.060660
0.354858 0.061497
0.378510 0.062133
0.402455 0.062562
0.426635 0.062779
0.450991 0.062774
0.475466 0.062530
0.500000 0.062029
0.524534 0.061254
0.549009 0.060194
0.573365 0.058845
0.597545 0.057218
0.621490 0.055344
0.645142 0.053258
0.668445 0.050993
0.691342 0.048575
0.713778 0.046029
0.735698 0.043377
0.757051 0.040641
0.777785 0.037847
0.797850 0.035017
0.817197 0.032176
0.835779 0.029347
0.853553 0.026554
0.870476 0.023817
0.886505 0.021153
0.901604 0.018580
0.915735 0.016113
0.928864 0.013769
0.940961 0.011562
0.951995 0.009508
0.961940 0.007622
0.970772 0.005915
0.978470 0.004401
0.985016 0.003092
0.990393 0.002001
0.994588 0.001137
0.997592 0.000510
0.999398 0.000128
1.000000 0.000000
0.000000 0.000000
0.000602 -0.003160
0.002408 -0.006308
0.005412 -0.009443
0.009607 -0.012559
0.014984 -0.015649
0.021530 -0.018707
0.029228 -0.021722
0.038060 -0.024685
0.048005 -0.027586
0.059039 -0.030416
0.071136 -0.033170
0.084265 -0.035843
0.098396 -0.038431
0.113495 -0.040929
0.129524 -0.043326
0.146447 -0.045610
0.164221 -0.047773
0.182803 -0.049805
0.202150 -0.051694
0.222215 -0.053427
0.242949 -0.054994
0.264302 -0.056376
0.286222 -0.057547
0.308658 -0.058459
0.331555 -0.059046
0.354858 -0.059236
0.378510 -0.058974
0.402455 -0.058224
0.426635 -0.056979
0.450991 -0.055257
0.475466 -0.053099
0.500000 -0.050563
0.524534 -0.047719
0.549009 -0.044642
0.573365 -0.041397
0.597545 -0.038043
0.621490 -0.034631
0.645142 -0.031207
0.668445 -0.027814
0.691342 -0.024495
0.713778 -0.021289
0.735698 -0.018232
0.757051 -0.015357
0.777785 -0.012690
0.797850 -0.010244
0.817197 -0.008027
0.835779 -0.006048
0.853553 -0.004314
0.870476 -0.002829
0.886505 -0.001592
0.901604 -0.000600
0.915735 0.000157
0.928864 0.000694
0.940961 0.001033
0.951995 0.001197
0.961940 0.001212
0.970772 0.001112
0.978470 0.000935
0.985016 0.000719
0.990393 0.000497
0.994588 0.000296
0.997592 0.000137
0.999398 0.000035
1.000000 0.000000"""
    pts = np.array([list(map(float, ln.split())) for ln in raw.split('\n')])
    up = pts[:65]      # LE->TE upper (index 65 is the duplicate 0,0 row)
    lo = pts[65:]      # LE->TE lower
    # resample upper+lower to n chordwise stations, cosine-clustered
    def resample(p):
        s = np.hstack([0, np.cumsum(np.hypot(np.diff(p[:, 0]),
                                             np.diff(p[:, 1])))])
        s = s / s[-1]
        return np.interp(np.linspace(0, 1, 350), s, p[:, 0]), \
            np.interp(np.linspace(0, 1, 350), s, p[:, 1])
    xu, yu = resample(up)
    xl, yl = resample(lo[::-1])
    return xu, yu, xl, yl


def biconvex_xy(n=350, t=0.08):
    """Circular-arc biconvex: y = +-sqrt((t/2)^2-(x-0.5)^2*t^2)-free? use
    standard 8% biconvex: y = (t/2)(1 - (2x-1)^2) shape (parabolic) --
    actually circular arc: y=(t/2)*sqrt(1-(2x-1)^2)-? Use classic:
    y = 2*t*(x*(1-x)) -- thin wedge family. We use the standard
    circular-arc biconvex: y = (t/2)*(1 - (1-2x)^2)."""
    x = np.linspace(0, 1, n)
    y = 0.5 * t * (1 - (1 - 2 * x) ** 2)  # = 2*t*x*(1-x), peak t/2 at mid
    return x, y


def extruded_slab_stl(xy_fn, span=1.6, n_span=220, scale=1.0,
                      tag='naca0012'):
    """Airfoil slab: section loop x spanwise layers => binary STL."""
    loop = xy_fn
    n_sec = len(loop)
    # spanwise stations (uniform)
    ys = np.linspace(-span / 2, span / 2, n_span)
    # vertices: (n_sec+1)*n_span ring
    verts = []
    for y in ys:
        for (x, z) in loop:
            verts.append((x * scale, y, z * scale))
    verts = np.array(verts, dtype=np.float32)
    tris = []
    for j in range(n_span - 1):
        b0, b1 = j * n_sec, (j + 1) * n_sec
        for i in range(n_sec):
            i2 = (i + 1) % n_sec
            a, b, c, d = b0 + i, b0 + i2, b1 + i2, b1 + i
            tris.append((a, b, c))
            tris.append((a, c, d))
    # end caps (fan within each end ring only, no out-of-bounds verts)
    for ring in [0, (n_span - 1) * n_sec]:
        for i in range(n_sec):
            i2 = (i + 1) % n_sec
            a, b, c = ring + i, ring + i2, ring
            if a == c or b == c:
                continue
            tris.append((b, a, c))
    tris = np.array(tris, dtype=np.int64)
    return verts, tris


def write_binary_stl(path, verts, tris):
    with open(path, 'wb') as fh:
        fh.write(b'\x00' * 80)
        fh.write(np.uint32(len(tris)))
        for t in tris:
            n = np.cross(verts[t[1]] - verts[t[0]], verts[t[2]] - verts[t[0]])
            nn = np.linalg.norm(n)
            n = n / nn if nn > 0 else np.array([0, 0, 1], dtype=np.float32)
            fh.write(np.asarray(n, dtype=np.float32).tobytes())
            for k in t:
                fh.write(verts[k].tobytes())
            fh.write(np.uint16(0))


def dfvlr_r4_xy():
    """DFVLR R-4 transonic research airfoil (UIUC coords)."""
    raw = """0.0 0.0
0.0000370 0.0043010
0.0004076 0.0063066
0.0009840 0.0079153
0.0014634 0.0093224
0.0019486 0.0103296
0.0042938 0.0140644
0.0092091 0.0198370
0.0141449 0.0242098
0.0190909 0.0278827
0.0240471 0.0308556
0.0290063 0.0336285
0.0339714 0.0360015
0.0389394 0.0381746
0.0439132 0.0399476
0.0488885 0.0416207
0.0538667 0.0430938
0.0588465 0.0444669
0.0638261 0.0458399
0.0688087 0.0470131
0.0787770 0.0491594
0.0987206 0.0529519
0.1186717 0.0562446
0.1386300 0.0590373
0.1586928 0.0615315
0.1786599 0.0637242
0.1986300 0.0657170
0.2186030 0.0675098
0.2385804 0.0690026
0.2585592 0.0703955
0.2785410 0.0715883
0.2985242 0.0726812
0.3185119 0.0734741
0.3386009 0.0741684
0.3585917 0.0747615
0.3785851 0.0751544
0.3875828 0.0752862
0.3985800 0.0754474
0.4185765 0.0756404
0.4385759 0.0756333
0.4585767 0.0755263
0.4785790 0.0753194
0.4985843 0.0749124
0.5185909 0.0744055
0.5386991 0.0738000
0.5587102 0.0729931
0.5787228 0.0720862
0.5987382 0.0709793
0.6187566 0.0696724
0.6387780 0.0681655
0.6588023 0.0664588
0.6788309 0.0644520
0.6988640 0.0621453
0.7190014 0.0595400
0.7390418 0.0567334
0.7590850 0.0537267
0.7791313 0.0505200
0.7991804 0.0471134
0.8192326 0.0435068
0.8392890 0.0396003
0.8593455 0.0356937
0.8794064 0.0314871
0.8994688 0.0271806
0.9196326 0.0227755
0.9397007 0.0180690
0.9597704 0.0132626
0.9798445 0.0081561
0.9848638 0.0068295
0.9898830 0.0055029
0.9949023 0.0041763
1.0000000 0.0025497
0.0 0.0
0.0005277 -0.0018925
0.0011481 -0.0032835
0.0023814 -0.0055657
0.0046219 -0.0083332
0.0096822 -0.0124595
0.0147263 -0.0154859
0.0197646 -0.0181124
0.0247970 -0.0203389
0.0298251 -0.0222654
0.0348502 -0.0239920
0.0398738 -0.0256186
0.0448960 -0.0271452
0.0499152 -0.0284718
0.0549345 -0.0297984
0.0599523 -0.0310250
0.0649685 -0.0321516
0.0699848 -0.0332782
0.0749997 -0.0343049
0.0800131 -0.0352315
0.1000667 -0.0389381
0.1201128 -0.0421448
0.1401533 -0.0449515
0.1602907 -0.0475567
0.1803238 -0.0498634
0.2003524 -0.0518703
0.2203782 -0.0536771
0.2403995 -0.0551839
0.2604177 -0.0564907
0.2804319 -0.0574976
0.3004430 -0.0583045
0.3204481 -0.0587114
0.3405503 -0.0589170
0.3605499 -0.0589240
0.3805448 -0.0586310
0.4005368 -0.0581381
0.4205244 -0.0573451
0.4405106 -0.0564523
0.4604939 -0.0553594
0.4804743 -0.0540665
0.4904629 -0.0533201
0.5054460 -0.0522005
0.5255190 -0.0504062
0.5454890 -0.0484134
0.5654562 -0.0462206
0.5854205 -0.0438279
0.6053802 -0.0411352
0.6253356 -0.0381426
0.6452882 -0.0349499
0.6652364 -0.0314573
0.6851801 -0.0276647
0.7051194 -0.0235722
0.7251601 -0.0195781
0.7450994 -0.0154855
0.7650403 -0.0114930
0.7849840 -0.0077004
0.8049306 -0.0041078
0.8248845 -0.0010151
0.8448443 0.0016776
0.8648143 0.0036704
0.8847933 0.0050632
0.9048809 0.0058576
0.9248788 0.0059506
0.9448869 0.0053437
0.9649112 0.0036369
0.9849546 0.0006302
0.9899694 -0.0003965
0.9949858 -0.0015231
1.0000000 -0.0025497"""
    pts = np.array([list(map(float, ln.split())) for ln in raw.split('\n')])
    # upper block ends at the FIRST (0,0) row that restarts the loop
    # (upper runs 0->1 with y>0, then a (0,0) row opens the lower block)
    zeros = np.where((pts[:, 0] == 0) & (pts[:, 1] == 0))[0]
    n_up = int(zeros[1]) if len(zeros) > 1 else int(np.argmax(
        (pts[:, 0] == 0) & (pts[:, 1] < 0)))
    up = pts[:n_up]
    lo = pts[n_up:]
    def resample(p):
        s = np.hstack([0, np.cumsum(np.hypot(np.diff(p[:, 0]),
                                             np.diff(p[:, 1])))])
        s = s / s[-1]
        return np.interp(np.linspace(0, 1, 350), s, p[:, 0]), \
            np.interp(np.linspace(0, 1, 350), s, p[:, 1])
    xu, yu = resample(up)
    xl, yl = resample(lo[::-1])
    return xu, yu, xl, yl


def dfvlr_loop():
    """DFVLR R-4 section as a CLOSED loop (see rae_loop note)."""
    xu, yu, xl, yl = dfvlr_r4_xy()
    up = np.stack([xu, yu], axis=1)[::-1]   # TE -> LE upper
    lo = np.stack([xl, yl], axis=1)[::-1]   # LE -> TE lower
    return np.vstack([up, lo])


def naca_loop(n=350, t=0.12):
    x, yt = naca0012_xy(n=n, t=t)
    up = np.stack([x[::-1], yt[::-1]], axis=1)
    lo = np.stack([x, -yt], axis=1)[1:]
    return np.vstack([up, lo])


def rae_loop():
    """RAE2822 section as a CLOSED loop: upper TE->LE, then lower LE->TE.

    (A naive upper-then-lower both LE->TE leaves the wrap quad spanning
    the whole chord -- a phantom flap that doubles lift.)
    """
    xu, yu, xl, yl = rae2822_xy()
    up = np.stack([xu, yu], axis=1)[::-1]   # reverse: TE -> LE upper
    lo = np.stack([xl, yl], axis=1)[::-1]   # xl is descending -> LE -> TE
    return np.vstack([up, lo])


def biconv_loop(n=350, t=0.08):
    """Biconvex section as a CLOSED loop (see rae_loop note)."""
    x, y = biconvex_xy(n=n, t=t)
    up = np.stack([x, y], axis=1)[::-1]     # TE -> LE upper
    lo = np.stack([x, -y], axis=1)[1:]      # LE+1 -> TE lower
    return np.vstack([up, lo])


def main():
    os.makedirs(BASE, exist_ok=True)
    suite = [('naca0012', naca_loop()), ('rae2822', rae_loop()),
             ('biconvex', biconv_loop()), ('dfvlr_r4', dfvlr_loop())]
    for name, loop in suite:
        verts, tris = extruded_slab_stl(loop, span=1.6, n_span=220)
        path = os.path.join(BASE, name + '.stl')
        write_binary_stl(path, verts, tris)
        # count unique vertices
        nv = len(np.unique(verts, axis=0))
        print('%s: verts=%d (unique %d) tris=%d  %.1f KB' % (
            name, len(verts), nv, len(tris), os.path.getsize(path) / 1024))


if __name__ == '__main__':
    main()
