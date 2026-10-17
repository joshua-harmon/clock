# Where the Pleiades sit in the Orion (Local) Arm: face-on around the Sun, and a cross-section through the arm.
# Arms: logarithmic spirals after Reid et al. (2019, ApJ 885, 131), as in make_orion_arm_map.py; the bands are the arm centre
# lines plus and minus the fitted arm width. Pleiades: l = 166.6, b = -23.5, 136 pc (444 ly; Lodieu et al. 2019, Gaia DR2).
# Sun's height above the galactic mid-plane: 20.8 pc = 68 ly (Bennett & Bovy 2019).
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
here = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(here, '..', '..', 'time', 'u6', 'u6_docs', 'pleiades')
KPC = 3261.56; R0 = 8.15
ARMS = {'Sagittarius–Carina Arm': (24, 6.04, 17.1, 1.0, 0.27, '#8fb4ff'),
        'Local (Orion) Arm':      (9, 8.26, 11.4, 11.4, 0.31, '#ffd36e'),
        'Perseus Arm':            (40, 8.87, 10.3, 8.7, 0.35, '#b8d0ff')}
def Rarm(beta, bk, rk, p1, p2):
    psi = p1 if beta < bk else p2
    return rk*math.exp(-math.radians(beta - bk)*math.tan(math.radians(psi)))
def arm(bk, rk, p1, p2, off=0.0, b0=-40, b1=40, n=800):
    beta = np.linspace(b0, b1, n); R = np.array([Rarm(b, bk, rk, p1, p2) for b in beta]) + off
    return R*np.sin(np.radians(beta))*KPC, (R*np.cos(np.radians(beta)) - R0)*KPC
def pos(l, b, d):
    lr, br = math.radians(l), math.radians(b); dp = d*math.cos(br)
    return dp*math.sin(lr), -dp*math.cos(lr), d*math.sin(br)
def galR(x, y):
    gx, gy = x/KPC, y/KPC + R0
    return math.hypot(gx, gy), math.degrees(math.atan2(gx, gy))
ZSUN = 68.0
P = pos(166.6, -23.5, 444)
LAND = [('Hyades', 180.1, -22.3, 153, 'cluster'), ('Orion Nebula', 209.0, -19.4, 1344, 'nebula'),
        ('Antares', 351.9, 15.1, 550, 'star'), ('Deneb', 84.3, 2.0, 2600, 'star'), ('Vela supernova remnant', 263.9, -3.3, 800, 'nebula'),
        ('Rosette Nebula', 206.3, -2.1, 5200, 'nebula'), ('Double Cluster', 134.6, -3.7, 7500, 'cluster')]
STY = {'star': ('o', '#ffcf6e', 6), 'cluster': ('D', '#9fe0ff', 6), 'nebula': ('s', '#ff8fb0', 5.5)}
la = ARMS['Local (Orion) Arm']
RP, bP = galR(P[0], P[1]); RS, bS = galR(0, 0)
dP = (RP - Rarm(bP, *la[:4]))*KPC; dS = (RS - Rarm(bS, *la[:4]))*KPC; W = la[4]*KPC
print('Pleiades %.0f ly from the arm centre line (negative = toward the galactic centre); Sun %.0f ly; arm width +/- %.0f ly' % (dP, dS, W))

# --- where they were 20 million years ago: straight-line trace-back from today's measured motion (relative to the Sun) ---
TBACK = 20.0                        # million years
T_ICRS_GAL = np.array([[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]])
KMS_LY_MYR = 3.33564                 # 1 km/s = 3.336 light years per million years
# name, RA, Dec (deg), parallax (mas), pmRA*cosDec, pmDec (mas/yr), radial velocity (km/s), colour, age note
MOV = [('Pleiades', 56.75, 24.12, 7.328, 20.01, -45.49, 5.8, '#ff7ad9', ''),        # medians of the Lodieu et al. 2019 members
       ('Hyades', 66.73, 15.87, 21.05, 101.0, -28.5, 39.4, '#9fe0ff', '')]            # Gaia Collaboration (Babusiaux et al.) 2018
# The Sun's own motion relative to the Local Standard of Rest (U toward the galactic centre, V toward l = 90, W north): Schönrich, Binney & Dehnen 2010
VSUN = np.array([11.1, 12.24, 7.25])
def helio(ra, de, plx, pmra, pmde, rv):
    a, d = math.radians(ra), math.radians(de); dist = 3261.56/plx
    r = np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
    e = np.array([-math.sin(a), math.cos(a), 0.0]); n = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
    v = rv*r + 4.74047*pmra/plx*e + 4.74047*pmde/plx*n            # km/s, ICRS
    pg, vg = T_ICRS_GAL @ (dist*r), T_ICRS_GAL @ v               # galactic X (to l=0), Y (to l=90), Z (north)
    return pg, vg
def mapxy(g): return g[1], -g[0]                                  # this map: x toward l = 90, y away from the galactic centre
PAST = []; VL = {}
for name, ra, de, plx, pmra, pmde, rv, col, note in MOV:
    pg, vg = helio(ra, de, plx, pmra, pmde, rv); vl = vg + VSUN                  # velocity relative to the Local Standard of Rest
    old = pg - TBACK*KMS_LY_MYR*vl
    VL[name] = vl; PAST.append((name, mapxy(pg), mapxy(old), old[2], pg[2], col, np.linalg.norm(vl), note))
    print('%s: %.1f km/s through the LSR; 20 Myr ago at (%.0f, %.0f) ly on the map' % ((name, np.linalg.norm(vl)) + mapxy(old)))
# Ursa Major moving group core (Merak, Phecda, Megrez, Alioth, Mizar, Alcor, 78 UMa): SIMBAD values (Hipparcos van Leeuwen 2007 / Gaia DR3), averaged
UMA = []
for ln in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ursa_major', 'uma_core_simbad.tsv'), encoding='utf-8').read().splitlines()[1:]:
    c = ln.split('\t'); UMA.append(helio(*[float(v) for v in c[1:7]]))
pg = np.mean([u[0] for u in UMA], axis=0); vg = np.mean([u[1] for u in UMA], axis=0); vl = vg + VSUN; old = pg - TBACK*KMS_LY_MYR*vl
UMACOL = '#8dffa8'
VL['Ursa Major core'] = vl; PAST.append(('Ursa Major core', mapxy(pg), mapxy(old), old[2], pg[2], UMACOL, np.linalg.norm(vl), ''))
print('Ursa Major core: %.0f ly away, %.1f km/s through the LSR; now (%.0f, %.0f, z %.0f); 20 Myr ago at (%.0f, %.0f, z %.0f)' % ((np.linalg.norm(pg), np.linalg.norm(vl)) + mapxy(pg) + (pg[2],) + mapxy(old) + (old[2],)))
SUNOLD = mapxy(-TBACK*KMS_LY_MYR*VSUN)
print('Sun: %.1f km/s through the LSR; 20 Myr ago at (%.0f, %.0f) ly' % ((np.linalg.norm(VSUN),) + SUNOLD))
# ---- curved paths instead of straight lines ----
# In the plane: epicycle (Hill) equations about the circular orbit at the Sun's distance, Oort constants A = 15.3, B = -11.9 km/s/kpc (Bovy 2017).
# Up and down: harmonic bounce through the mid-plane for a local density of 0.084 solar masses per cubic parsec (McKee, Parravano & Hollenbach 2015).
KMS = 1.0227; PC = 3.26156                                   # 1 km/s = 1.0227 pc per million years; 1 pc = 3.26156 ly
OORT_A, OORT_B = 15.3, -11.9
OM = (OORT_A - OORT_B)*KMS/1000; AA = OORT_A*KMS/1000        # per million years
RHO0 = 0.084; NU = math.sqrt(4*math.pi*4.4985e-3*RHO0)      # G = 4.4985e-3 pc^3 / (solar mass Myr^2)
print('epicycle period %.0f Myr, vertical period %.0f Myr' % (2*math.pi/math.sqrt(-4*OM*OORT_B*KMS/1000), 2*math.pi/NU))
def orbit(now_xy, zmid, v, T=TBACK, n=400):
    """now_xy: map position (ly); zmid: height above mid-plane (ly); v: LSR velocity (U to GC, V rotation, W north, km/s).
    Returns map x, map y, height (ly) from now (index 0) back to T million years ago (index n)."""
    st = np.array([now_xy[1]/PC, now_xy[0]/PC, -v[0]*KMS, v[1]*KMS])          # x outward, y along rotation (pc, pc/Myr)
    def fd(q): return np.array([q[2], q[3], 2*OM*q[3] + 4*OM*AA*q[0], -2*OM*q[2]])
    dt = -T/n; out_ = [st.copy()]
    for _ in range(n):
        k1 = fd(st); k2 = fd(st + dt/2*k1); k3 = fd(st + dt/2*k2); k4 = fd(st + dt*k3); st = st + dt/6*(k1 + 2*k2 + 2*k3 + k4); out_.append(st.copy())
    Q = np.array(out_); t = np.linspace(0, -T, n + 1)
    z = (zmid/PC)*np.cos(NU*t) + (v[2]*KMS/NU)*np.sin(NU*t)
    return Q[:, 1]*PC, Q[:, 0]*PC, z*PC
VL['Sun'] = VSUN.copy(); PATH = {}; _new = []
for name, now, old, oz, nz, col, spd, note in PAST:
    PATH[name] = orbit(now, nz + ZSUN, VL[name])
    px, py, pz = PATH[name]; _new.append((name, now, (px[-1], py[-1]), pz[-1] - ZSUN, nz, col, spd, note))
    print('%s curved: 20 Myr ago at (%.0f, %.0f), height %.0f (straight line gave (%.0f, %.0f), %.0f)' % (name, px[-1], py[-1], pz[-1], old[0], old[1], oz + ZSUN))
PAST = _new
PATH['Sun'] = orbit((0.0, 0.0), ZSUN, VSUN)
print('Sun curved: 20 Myr ago at (%.0f, %.0f), height %.0f (straight line gave (%.0f, %.0f), %.0f)' % (PATH['Sun'][0][-1], PATH['Sun'][1][-1], PATH['Sun'][2][-1], SUNOLD[0], SUNOLD[1], ZSUN - TBACK*KMS_LY_MYR*VSUN[2]))
SUNOLD = (PATH['Sun'][0][-1], PATH['Sun'][1][-1])
def curve(axh, xs, ys, col, lw, ms, rx, ry, z=9, alpha=.9):
    """dashed path from 20 Myr ago to now, with an arrowhead pointing at the present position"""
    xs, ys = np.asarray(xs), np.asarray(ys)
    axh.plot(xs, ys, color=col, lw=lw, ls=(0, (4, 3)), alpha=alpha, zorder=z, solid_capstyle='round')
    d = np.hypot((xs - xs[0])/rx, (ys - ys[0])/ry)
    i1 = int(np.argmax(d >= 0.012)); i2 = int(np.argmax(d >= 0.03))
    if i2 > i1 > 0:
        axh.annotate('', xy=(xs[i1], ys[i1]), xytext=(xs[i2], ys[i2]), arrowprops=dict(arrowstyle='-|>', mutation_scale=ms, color=col, lw=lw, shrinkA=0, shrinkB=0), zorder=z)
BG = '#05070d'
fig = plt.figure(figsize=(20, 10.5), dpi=125); fig.patch.set_facecolor(BG)
ax = fig.add_axes([0.01, 0.11, 0.56, 0.80]); ax.set_facecolor(BG)
M = 3000
for name, (bk, rk, p1, p2, w, col) in ARMS.items():
    xi, yi = arm(bk, rk, p1, p2, -w); xo, yo = arm(bk, rk, p1, p2, +w); xc, yc = arm(bk, rk, p1, p2)
    ax.fill(np.r_[xi, xo[::-1]], np.r_[yi, yo[::-1]], color=col, alpha=.13 if 'Local' in name else .08, lw=0, zorder=1)
    ax.plot(xc, yc, color=col, lw=1.2, ls=(0, (6, 4)), alpha=.8, zorder=2)
    ax.plot(xi, yi, color=col, lw=.6, alpha=.35, zorder=2); ax.plot(xo, yo, color=col, lw=.6, alpha=.35, zorder=2)
    k = np.argmin(abs(xc + 0.70*M) + 1e6*(abs(yc) > M))
    if abs(yc[k]) < 0.9*M and abs(xc[k]) < M: ax.text(xc[k], yc[k] + 120, name + ' (centre line)', color=col, fontsize=12, style='italic', ha='center', zorder=6,
                               bbox=dict(boxstyle='round,pad=.2', fc=BG, ec='none', alpha=.75))
xi, yi = arm(*la[:4], -la[4]); xo, yo = arm(*la[:4], la[4])
ki = np.argmin(abs(xi + 0.45*M)); ko = np.argmin(abs(xo - 0.55*M))
ax.text(xi[ki], yi[ki] - 70, 'inner edge (toward the galactic centre)', color='#ffd36e', fontsize=9.5, ha='center', va='top', alpha=.9, zorder=6)
ax.text(xo[ko], yo[ko] + 70, 'outer edge', color='#ffd36e', fontsize=9.5, ha='center', va='bottom', alpha=.9, zorder=6)
for rr in (500, 1000, 2000, 3000):
    ax.add_patch(Circle((0, 0), rr, fill=False, ec='#d9c27a', lw=.5, ls=(0, (3, 4)), alpha=.45, zorder=3))
    ax.text(rr*0.71 + 30, -rr*0.71 - 30, '%s ly' % format(rr, ','), color='#cdb874', fontsize=8.5, zorder=7)
for l in range(0, 360, 30):
    lr = math.radians(l); ax.plot([0, 1.5*M*math.sin(lr)], [0, -1.5*M*math.cos(lr)], color='#d9c27a', lw=.4, alpha=.3, zorder=3)
    s, c = math.sin(lr), -math.cos(lr); t = min(M*.95/abs(s) if abs(s) > 1e-9 else 1e9, M*.95/abs(c) if abs(c) > 1e-9 else 1e9)
    ax.text(t*s, t*c, '%d°' % l, color='#f0d98a', fontsize=9, ha='center', va='center', bbox=dict(boxstyle='round,pad=.15', fc=BG, ec='none'), zorder=7)
for name, l, b, d, kind in LAND:
    x, y, z = pos(l, b, d); mk, col, ms = STY[kind]
    if abs(x) < M and abs(y) < M:
        ax.plot(x, y, marker=mk, ms=ms, color=col, ls='none', zorder=8)
        ax.annotate(name, (x, y), xytext=(-6 if name in ('Vela supernova remnant', 'Orion Nebula') else (14 if name=='Hyades' else 6), 6 if name in ('Hyades',) else 4), ha='right' if name in ('Vela supernova remnant', 'Orion Nebula') else 'left', textcoords='offset points', color=col, fontsize=9.5, zorder=9)
LBL = {'Pleiades': (1500, 1250), 'Alpha Persei cluster': (1750, 600), 'Hyades': (900, -2050)}
for name, now, old, oz, nz, col, spd, note in PAST:
    ax.add_patch(Circle(old, 38, fill=False, ec=col, lw=2.0, ls=(0, (3, 2)) if note else '-', zorder=10))
    curve(ax, PATH[name][0], PATH[name][1], col, 1.5, 20, 2*M, 2*M)
ax.add_patch(Circle(SUNOLD, 40, fill=False, ec='#ffe24a', lw=2.4, zorder=12))
curve(ax, PATH['Sun'][0], PATH['Sun'][1], '#ffe24a', 2.0, 22, 2*M, 2*M, z=12)
ax.plot(P[0], P[1], marker='o', ms=7, color='#ff7ad9', ls='none', zorder=11)
_u = [p for p in PAST if p[0] == 'Ursa Major core'][0]
ax.plot(*_u[1], marker='o', ms=6, color=UMACOL, ls='none', zorder=11)
_pp = PATH['Ursa Major core']; _m = (_pp[0][len(_pp[0])//2], _pp[1][len(_pp[1])//2])
ax.text(_m[0] + 70, _m[1], 'Ursa Major core\n(Big Dipper stars)', color=UMACOL, fontsize=9.5, ha='left', va='center', zorder=11)
ax.annotate('Pleiades', (P[0], P[1]), xytext=(0, 8), textcoords='offset points', color='#ff9be3', fontsize=9.5, ha='center', va='bottom', zorder=11)
ax.plot(0, 0, marker='o', ms=7, color='#ffe24a', ls='none', zorder=10); ax.annotate('Sun', (0, 0), xytext=(8, -16), textcoords='offset points', color='#ffe24a', fontsize=12, weight='bold')
ax.text(0, -M*.88, '↓ toward the galactic centre (26,600 ly)', color='#ddd', fontsize=11, ha='center')
ax.text(0, M*.86, 'away from the galactic centre ↑', color='#ddd', fontsize=11, ha='center')
ax.set_xlim(-M, M); ax.set_ylim(-M, M); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('From above (north galactic pole): 3,000 light years around the Sun, now and 20 million years ago', color='white', fontsize=13.5)

# cross-section through the Local Arm along the line from the galactic centre through the Pleiades
XOFF = {}
bx = fig.add_axes([0.60, 0.22, 0.38, 0.52]); bx.set_facecolor('#0a0f1c')
for sp in bx.spines.values(): sp.set_color('#3a4a60')
bx.tick_params(colors='#9aa', labelsize=8.5, length=3)
H = 600
bx.axvspan(-W, W, color='#ffd36e', alpha=.12, lw=0)
bx.axvline(0, color='#ffd36e', ls=(0, (6, 4)), lw=1.1); bx.axvline(-W, color='#ffd36e', lw=.6, alpha=.5); bx.axvline(W, color='#ffd36e', lw=.6, alpha=.5)
bx.axhline(0, color='#9fb6d8', lw=.8, ls=(0, (2, 3)))
# (text removed) bx.text(-W + 30, H*0.86, 'inner edge', color='#ffd36e', fontsize=9.5); bx.text(W - 30, H*0.86, 'outer edge', color='#ffd36e', fontsize=9.5, ha='right')
# (text removed) bx.text(15, H*0.72, 'arm centre line', color='#ffd36e', fontsize=9.5)
# (text removed) bx.text(1.55*W, 12, 'galactic mid-plane', color='#9fb6d8', fontsize=9, ha='right', va='bottom')
bx.plot(dP, P[2] + ZSUN, marker='o', ms=7, color='#ff7ad9', ls='none', zorder=7)
bx.annotate('Pleiades', (dP, P[2] + ZSUN), xytext=(10, -12), textcoords='offset points', color='#ff9be3', fontsize=9)

bx.plot(dS, ZSUN, marker='o', ms=7, color='#ffe24a', ls='none', zorder=7); bx.annotate('Sun', (dS, ZSUN), xytext=(-9, -4), textcoords='offset points', color='#ffe24a', fontsize=9, ha='right')
def armoff(xy):
    Rg, bg = galR(*xy); return (Rg - Rarm(bg, *la[:4]))*KPC
SUNOLDZ = PATH['Sun'][2][-1]
XPATH = {n: (np.array([armoff((a, b)) for a, b in zip(p[0], p[1])]), p[2]) for n, p in PATH.items()}
XS = [(armoff(SUNOLD), SUNOLDZ, dS, ZSUN, '#ffe24a', 'Sun')]
for name, now, old, oz, nz, col, spd, note in PAST:
    XS.append((armoff(old), oz + ZSUN, armoff(now), nz + ZSUN, col, name))
for ox, oy, nx, ny, col, name in XS:
    print('cross-section %s: then (%.0f, %.0f) now (%.0f, %.0f)' % (name, ox, oy, nx, ny))
    bx.scatter([ox], [oy], s=65 if name == 'Sun' else 50, facecolors='none', edgecolors=col, linewidths=1.8 if name == 'Sun' else 1.5, zorder=6)
    curve(bx, XPATH[name][0], XPATH[name][1], col, 1.4, 18, 5000, 1200, z=5)
    if name in ('Hyades', 'Ursa Major core'):
        bx.plot(nx, ny, marker='o', ms=6, color=col, ls='none', zorder=6)
        bx.annotate(name, (nx, ny), xytext=(-6, -30) if name == 'Hyades' else (-9, 2), ha='right', textcoords='offset points', color=col, fontsize=9)
allx = [v for p in XPATH.values() for v in p[0]] + [-1.6*W, 1.6*W]; ally = [v for p in XPATH.values() for v in p[1]] + [-H, H]
bx.set_xlim(min(allx) - 250, max(allx) + 250); bx.set_ylim(min(ally) - 120, max(ally) + 120)
bx.set_xlabel('← galactic centre        across the arm (ly)        outward →', color='#ccc', fontsize=10)
bx.set_ylabel('height above mid-plane (ly)', color='#ccc', fontsize=10)
# (text removed) bx.set_title('Cross-section through the Orion Arm at the Pleiades', color='white', fontsize=13)
fig.suptitle('The Pleiades in the Orion Arm', color='white', fontsize=17, y=0.975)
# (text removed) fig.text(0.60, 0.77, 'The Pleiades lie inside the Orion Arm, in its inner half: about %s ly from the arm’s centre line toward\nthe galactic centre, '
# (text removed)          'within the arm’s fitted width of ±%s ly. In this model the Sun lies a little farther in, just beyond the inner edge.' % (format(int(round(-dP, -1)), ','), format(int(round(W, -1)), ',')),
# (text removed)          color='#e8e8e8', fontsize=11, va='bottom')
fig.text(0.01, 0.005, 'Arms after Reid et al. 2019 (ApJ 885, 131; maser parallaxes): dashed = arm centre line, band = centre ± fitted width; arm edges are not sharp, and other tracers place the arm slightly differently. '
         'Pleiades distance 136 pc = 444 ly (Lodieu et al. 2019, Gaia DR2). Sun 68 ly above the mid-plane (Bennett & Bovy 2019).\n'
         'Circles: where the Sun and these stars and clusters were 20 million years ago, traced back along curved paths from their measured proper motion, parallax and radial velocity (Gaia DR2), in the frame that circles the galaxy at the Sun’s distance (the Local Standard of Rest; Sun’s own motion from Schönrich et al. 2010), drawn around where the Sun is now.\n'
         'Paths: epicycles in the plane (Oort constants A = 15.3, B = −11.9 km/s/kpc; Bovy 2017) and a bounce through the mid-plane every 91 million years (local density 0.084 solar masses per cubic parsec; McKee et al. 2015).\nThe band shows the arm as it is now. Ursa Major core = average of Merak, Phecda, Megrez, Alioth, Mizar, Alcor and 78 UMa (SIMBAD: Hipparcos, van Leeuwen 2007; Gaia DR3).',
         color='#9aa', fontsize=9.5, va='bottom')
fig.savefig(os.path.join(out, 'pleiades_in_orion_arm.png'), facecolor=BG, bbox_inches='tight'); print('saved')

# ---- second figure: two side views, height above the galactic mid-plane, now (dots) and 20 million years ago (circles) ----
SIDEOFF = {(0, 'Ursa Major core'): (-2, 14), (0, 'Pleiades'): (8, -14), (0, 'Hyades'): (-10, 4), (1, 'Pleiades'): (8, -14), (1, 'Sun'): (-8, 6)}
SIDE = [('Sun', (0.0, 0.0), SUNOLD, ZSUN, SUNOLDZ, '#ffe24a')] + [(n, now, old, nz + ZSUN, oz + ZSUN, col) for n, now, old, oz, nz, col, spd, note in PAST]
fig2 = plt.figure(figsize=(18, 7.2), dpi=125); fig2.patch.set_facecolor(BG)
VIEWS = [(1, '← away from the galactic centre        light years        toward it →', 'Seen edge-on, from behind, looking forward along the Sun’s orbit'),
         (0, '← behind        light years, along the direction of galactic rotation        ahead →', 'Seen edge-on, looking out from the galactic centre')]
for k, (ix, xl, tt) in enumerate(VIEWS):
    sx = fig2.add_axes([0.05 + 0.49*k, 0.13, 0.44, 0.74]); sx.set_facecolor('#0a0f1c')
    for sp in sx.spines.values(): sp.set_color('#3a4a60')
    sx.tick_params(colors='#9aa', labelsize=9)
    sx.axhline(0, color='#9fb6d8', lw=.9, ls=(0, (2, 3)))
    xs = []
    for name, now, old, hz, oz, col in SIDE:
        nx, ox = now[ix], old[ix]; xs += [nx, ox]
        curve(sx, PATH[name][ix], PATH[name][2], col, 1.4, 18, 4000, 1200, z=5); xs += list(PATH[name][ix])
        sx.scatter([ox], [oz], s=50, facecolors='none', edgecolors=col, linewidths=1.6, zorder=6)
        sx.plot(nx, hz, marker='o', ms=7, color=col, ls='none', zorder=7)
        sx.annotate(name, (nx, hz), xytext=SIDEOFF.get((k, name), (8, 6)), textcoords='offset points', color=col, fontsize=9.5,
                    ha='right' if SIDEOFF.get((k, name), (8, 6))[0] < 0 else 'left', zorder=8)
    sx.set_xlim(min(xs) - 300, max(xs) + 300); sx.set_ylim(-600, 600)
    if ix == 1: sx.invert_xaxis()
    sx.set_xlabel(xl, color='#ccc', fontsize=10); sx.set_ylabel('height above mid-plane (ly)', color='#ccc', fontsize=10)
    sx.set_title(tt, color='white', fontsize=12.5)
fig2.text(0.5, 0.955, 'Height above the galactic mid-plane: now (dots) and 20 million years ago (circles)', color='white', fontsize=15, ha='center')
fig2.text(0.01, 0.01, 'Same data and frame as pleiades_in_orion_arm.png (Local Standard of Rest; curved paths: epicycle in the plane, bounce through the mid-plane). Dotted line = galactic mid-plane. '
          'The vertical scale is stretched compared with the horizontal.', color='#9aa', fontsize=9.5)
fig2.savefig(os.path.join(out, 'pleiades_side_views.png'), facecolor=BG, bbox_inches='tight'); print('saved side views')

# ---- third figure: velocity diagram (relative to the Local Standard of Rest), every arrow from one common point ----
KV = TBACK*KMS_LY_MYR
VOFF = {(0, 'Pleiades'): (-8, -16), (1, 'Pleiades'): (-8, -16)}
VEL = {n: np.array([VL[n][1], -VL[n][0], VL[n][2]]) for n, now, old, hz, oz, col in SIDE}   # km/s: (rotation, away from GC, north)
VCOL = {n: col for n, now, old, hz, oz, col in SIDE}
for n, v in VEL.items(): print('velocity %s: %.1f %.1f %.1f  |v| %.1f' % (n, v[0], v[1], v[2], np.linalg.norm(v)))
def ang(a, b): return math.degrees(math.acos(np.dot(a, b)/np.linalg.norm(a)/np.linalg.norm(b)))
PAIRS = [('Sun', 'Ursa Major core'), ('Hyades', 'Pleiades')]
for a, b in PAIRS: print('%s vs %s: %.0f deg apart, %.1f km/s relative' % (a, b, ang(VEL[a], VEL[b]), np.linalg.norm(VEL[a] - VEL[b])))
fig3 = plt.figure(figsize=(16, 7.6), dpi=125); fig3.patch.set_facecolor(BG)
PAN = [((0, 1), '← behind     along galactic rotation (km/s)     ahead →', '← toward galactic centre   (km/s)   away →', 'From above (north galactic pole)'),
       ((0, 2), '← behind     along galactic rotation (km/s)     ahead →', '← south   (km/s)   north →', 'Edge-on, looking out from the galactic centre')]
for k, ((i, j), xl, yl, tt) in enumerate(PAN):
    vx = fig3.add_axes([0.05 + 0.49*k, 0.14, 0.42, 0.72]); vx.set_facecolor('#0a0f1c')
    for sp in vx.spines.values(): sp.set_color('#3a4a60')
    vx.tick_params(colors='#9aa', labelsize=9)
    vx.axhline(0, color='#3a4a60', lw=.8); vx.axvline(0, color='#3a4a60', lw=.8)
    for r in (10, 20, 30): vx.add_patch(Circle((0, 0), r, fill=False, ec='#3a4a60', lw=.6, ls=(0, (2, 3))))
    for a, b in PAIRS:
        vx.plot([VEL[a][i], VEL[b][i]], [VEL[a][j], VEL[b][j]], color='#bbbbbb', lw=.9, ls=(0, (1, 2)), zorder=3)
    for n, v in VEL.items():
        vx.annotate('', xy=(v[i], v[j]), xytext=(0, 0), arrowprops=dict(arrowstyle='-|>', mutation_scale=20, color=VCOL[n], lw=2.0, shrinkA=0, shrinkB=0), zorder=5)
        off = VOFF.get((k, n), (6, 4))
        vx.annotate(n, (v[i], v[j]), xytext=off, textcoords='offset points', color=VCOL[n], fontsize=10, ha='right' if off[0] < 0 else 'left', zorder=6)
    vx.plot(0, 0, marker='o', ms=5, color='white', zorder=6)
    vx.set_xlim(-36, 36); vx.set_ylim(-36, 36); vx.set_aspect('equal')
    vx.set_xlabel(xl, color='#ccc', fontsize=10); vx.set_ylabel(yl, color='#ccc', fontsize=10); vx.set_title(tt, color='white', fontsize=12.5)
txt = '   '.join('%s and %s: %d° apart, %d km/s relative speed' % (a, b, round(ang(VEL[a], VEL[b])), round(np.linalg.norm(VEL[a] - VEL[b]))) for a, b in PAIRS)
fig3.text(0.5, 0.955, 'Velocities: every arrow drawn from the same starting point', color='white', fontsize=15, ha='center')
fig3.text(0.5, 0.905, txt, color='#e8e8e8', fontsize=11, ha='center')
fig3.text(0.01, 0.01, 'White dot = the Local Standard of Rest (a circular orbit around the galaxy at the Sun’s distance, about 220 km/s, which is subtracted). Rings at 10, 20, 30 km/s. '
          'Dotted lines join the pairs being compared. Same data as pleiades_in_orion_arm.png.', color='#9aa', fontsize=9.5)
fig3.savefig(os.path.join(out, 'pleiades_velocities.png'), facecolor=BG, bbox_inches='tight'); print('saved velocities')

# ---- fourth figure: the same curved paths followed back 120 million years (about the age of the Pleiades) ----
LOFF = {'Sun': (-70, 22), 'Ursa Major core': (-70, 48), 'Hyades': (-70, -32), 'Pleiades': (70, -40)}
TLONG = 120.0; NL = 2400
LONG = {'Sun': orbit((0.0, 0.0), ZSUN, VSUN, T=TLONG, n=NL)}
for name, now, old, oz, nz, col, spd, note in PAST: LONG[name] = orbit(now, nz + ZSUN, VL[name], T=TLONG, n=NL)
LCOL = {'Sun': '#ffe24a'}; LCOL.update({p[0]: p[5] for p in PAST})
ORDER = ['Sun', 'Pleiades', 'Hyades', 'Ursa Major core']
tl = np.linspace(0, -TLONG, NL + 1)
for n in ORDER: print('%s 120 Myr ago: (%.0f, %.0f) height %.0f; height range %.0f to %.0f' % (n, LONG[n][0][-1], LONG[n][1][-1], LONG[n][2][-1], LONG[n][2].min(), LONG[n][2].max()))
fig4 = plt.figure(figsize=(16, 13.5), dpi=125); fig4.patch.set_facecolor(BG)
# left: from above, relative to the circular orbit at the Sun's distance
gx = fig4.add_axes([0.06, 0.50, 0.90, 0.42]); gx.set_facecolor('#0a0f1c')
for sp in gx.spines.values(): sp.set_color('#3a4a60')
gx.tick_params(colors='#9aa', labelsize=9)
gx.axhline(0, color='#3a4a60', lw=.7); gx.axvline(0, color='#3a4a60', lw=.7)
allxy = []
for n in ORDER:
    px, py, pz = LONG[n]; col = LCOL[n]; allxy += [px, py]
    gx.plot(px, py, color=col, lw=1.6 if n == 'Sun' else 1.3, alpha=.9, zorder=5)
    for tk in range(20, int(TLONG) + 1, 20):
        i = int(round(tk/TLONG*NL)); gx.plot(px[i], py[i], marker='o', ms=3.5, color=col, zorder=6)
        if n == 'Sun': gx.annotate('%d' % tk, (px[i], py[i]), xytext=(5, 4), textcoords='offset points', color=col, fontsize=8, alpha=.85)
    gx.scatter([px[-1]], [py[-1]], s=50, facecolors='none', edgecolors=col, linewidths=1.6, zorder=7)
    gx.plot(px[0], py[0], marker='o', ms=7, color=col, ls='none', zorder=8)
    gx.annotate(n, (px[0], py[0]), xytext=LOFF.get(n, (8, 6)), textcoords='offset points', color=col, fontsize=10, ha='right' if LOFF.get(n, (8, 6))[0] < 0 else 'left', va='center', zorder=9,
                arrowprops=dict(arrowstyle='-', color=col, lw=.6, alpha=.8, shrinkA=2, shrinkB=4))
lim = max(abs(np.concatenate(allxy)).max()*1.08, 1000)
ay = np.concatenate([LONG[n][1] for n in ORDER]); gx.set_xlim(-lim, lim); gx.set_ylim(ay.min() - 1800, ay.max() + 1800); gx.set_aspect('equal')
gx.set_xlabel('← behind        along galactic rotation (ly)        ahead →', color='#ccc', fontsize=10)
gx.set_ylabel('← toward galactic centre        (ly)        away →', color='#ccc', fontsize=10)
gx.set_title('From above, riding along with the circular orbit at the Sun’s distance\n(small dots every 20 million years; numbers on the Sun’s path)', color='white', fontsize=11.5)
# right: height above the mid-plane against time
hx = fig4.add_axes([0.06, 0.07, 0.86, 0.34]); hx.set_facecolor('#0a0f1c')
for sp in hx.spines.values(): sp.set_color('#3a4a60')
hx.tick_params(colors='#9aa', labelsize=9)
hx.axhline(0, color='#9fb6d8', lw=.9, ls=(0, (2, 3)))
for tk in range(-120, 1, 20): hx.axvline(tk, color='#3a4a60', lw=.5, alpha=.6)
for n in ORDER:
    pz = LONG[n][2]; col = LCOL[n]
    hx.plot(tl, pz, color=col, lw=2.0 if n == 'Sun' else 1.5, alpha=.95)
    hx.plot(0, pz[0], marker='o', ms=7, color=col); hx.annotate(n, (0, pz[0]), xytext=(8, 0), textcoords='offset points', color=col, fontsize=10, va='center')
    hx.plot(-20, pz[int(round(20/TLONG*NL))], marker='o', ms=7, mfc='none', mec=col, mew=1.5)
hx.axvline(-20, color='#bbbbbb', lw=.8, ls=(0, (4, 3))); hx.text(-20, 1.0, ' 20 Myr ago', color='#bbb', fontsize=9, transform=hx.get_xaxis_transform(), va='top')
hx.set_xlim(-TLONG - 3, 14); hx.set_ylim(-500, 500)
hx.set_xlabel('millions of years from now', color='#ccc', fontsize=10); hx.set_ylabel('height above mid-plane (ly)', color='#ccc', fontsize=10)
hx.set_title('Up and down through the galactic mid-plane', color='white', fontsize=12.5)
fig4.text(0.5, 0.975, 'The last 120 million years (about the age of the Pleiades)', color='white', fontsize=15, ha='center')
fig4.text(0.01, 0.005, 'Same measurements and model as pleiades_in_orion_arm.png: epicycles in the plane (Oort constants A = 15.3, B = −11.9 km/s/kpc; Bovy 2017), a bounce through the mid-plane every 91 million years '
          '(local density 0.084 solar masses per cubic parsec; McKee et al. 2015).\nThis simple model gets less exact the further back it goes; real paths are also nudged by spiral arms and gas clouds. '
          'The spiral arms are not drawn, because over 120 million years the arm pattern moves a long way relative to the stars.', color='#9aa', fontsize=9.5)
fig4.savefig(os.path.join(out, 'pleiades_120_million_years.png'), facecolor=BG, bbox_inches='tight'); print('saved 120 Myr')

MOFF = {'Pleiades': 12, 'Hyades': 0, 'Ursa Major core': -13, 'Sun': -26}
MLAB = {'Sun': (-60, 40), 'Ursa Major core': (60, 45), 'Hyades': (-70, -30), 'Pleiades': (60, -45)}
# ---- fifth figure: motion around the galaxy held fixed; only in/out and up/down over 120 million years ----
fig5 = plt.figure(figsize=(16, 13.5), dpi=125); fig5.patch.set_facecolor(BG)
mx = fig5.add_axes([0.07, 0.50, 0.86, 0.42]); mx.set_facecolor('#0a0f1c')
tx = fig5.add_axes([0.07, 0.07, 0.86, 0.34]); tx.set_facecolor('#0a0f1c')
for a_ in (mx, tx):
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
    a_.tick_params(colors='#9aa', labelsize=9)
mx.axhline(0, color='#9fb6d8', lw=.9, ls=(0, (2, 3))); mx.axvline(0, color='#3a4a60', lw=.7)
tx.axhline(0, color='#3a4a60', lw=.7)
for tk in range(-120, 1, 20): tx.axvline(tk, color='#3a4a60', lw=.5, alpha=.6)
for n in ORDER:
    px, py, pz = LONG[n]; col = LCOL[n]
    mx.plot(py, pz, color=col, lw=1.8 if n == 'Sun' else 1.3, alpha=.9, zorder=5)
    for tk in range(20, int(TLONG) + 1, 20):
        i = int(round(tk/TLONG*NL)); mx.plot(py[i], pz[i], marker='o', ms=3.5, color=col, zorder=6)
        if n == 'Sun': mx.annotate('%d' % tk, (py[i], pz[i]), xytext=(5, 4), textcoords='offset points', color=col, fontsize=8.5)
    mx.scatter([py[-1]], [pz[-1]], s=50, facecolors='none', edgecolors=col, linewidths=1.6, zorder=7)
    mx.plot(py[0], pz[0], marker='o', ms=7, color=col, ls='none', zorder=8)
    tx.plot(tl, py, color=col, lw=2.0 if n == 'Sun' else 1.5, alpha=.95)
    tx.plot(0, py[0], marker='o', ms=7, color=col); tx.annotate(n, (0, py[0]), xytext=(8, MOFF.get(n, 0)), textcoords='offset points', color=col, fontsize=10, va='center')
for n, (dx, dy) in MLAB.items():
    mx.annotate(n, (LONG[n][1][0], LONG[n][2][0]), xytext=(dx, dy), textcoords='offset points', color=LCOL[n], fontsize=10, ha='right' if dx < 0 else 'left', va='center',
                arrowprops=dict(arrowstyle='-', color=LCOL[n], lw=.6, alpha=.8, shrinkA=2, shrinkB=4), zorder=9)
mx.set_ylim(-480, 480)
mx.set_xlabel('← toward the galactic centre        light years        away from it →', color='#ccc', fontsize=10)
mx.set_ylabel('height above mid-plane (ly)', color='#ccc', fontsize=10)
mx.set_title('Seen edge-on with the trip around the galaxy taken out: in-and-out against up-and-down\n(small dots every 20 million years; numbers on the Sun’s path; circles = 120 million years ago)', color='white', fontsize=11.5)
tx.set_xlim(-TLONG - 3, 22)
tx.set_xlabel('millions of years from now', color='#ccc', fontsize=10)
tx.set_ylabel('← toward galactic centre   (ly)   away →', color='#ccc', fontsize=10)
tx.set_title('In and out from the galactic centre over time', color='white', fontsize=12.5)
fig5.text(0.5, 0.975, 'The last 120 million years, holding the motion around the galaxy constant', color='white', fontsize=15, ha='center')
fig5.text(0.01, 0.005, 'Same paths as pleiades_120_million_years.png with the forward and backward drift along the orbit removed. Zero = the Sun’s present distance from the galactic centre and the galactic mid-plane. '
          '\nIn and out repeats about every 171 million years, up and down about every 91 million years, so the shapes do not close. The vertical scale is stretched compared with the horizontal.', color='#9aa', fontsize=9.5)
fig5.savefig(os.path.join(out, 'pleiades_120_million_years_in_out_up_down.png'), facecolor=BG, bbox_inches='tight'); print('saved meridional')

# ---- sixth figure: the same in-and-out / up-and-down view over a longer span ----
def place_labels(fg, items):
    """put each text near its point where it touches no drawn line and no other label"""
    fg.canvas.draw(); rd = fg.canvas.get_renderer(); boxes = []
    for axh, xy, txt, col, fs in items:
        segs = []
        for ln in axh.get_lines():
            d = ln.get_transform().transform(np.column_stack([np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)]))
            if len(d) > 1:
                for p0, p1 in zip(d[:-1], d[1:]):
                    k = max(2, int(np.hypot(*(p1 - p0))/2) + 1); segs.append(np.linspace(p0, p1, k))
            else: segs.append(d)
        pts = np.vstack(segs)
        best = None
        for r in (5, 9, 14, 20, 28, 38):
            for a in range(0, 360, 30):
                dx, dy = r*math.cos(math.radians(a)), r*math.sin(math.radians(a))
                t = axh.annotate(txt, xy, xytext=(dx, dy), textcoords='offset points', color=col, fontsize=fs, ha='left' if dx >= -1 else 'right', va='bottom' if dy >= -1 else 'top', zorder=9)
                bb = t.get_window_extent(rd).expanded(1.08, 1.15)
                hit = np.any((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)) or any(bb.overlaps(o) for o in boxes)
                inside = axh.get_window_extent(rd).contains(bb.x0, bb.y0) and axh.get_window_extent(rd).contains(bb.x1, bb.y1)
                if not hit and inside: best = t; break
                t.remove()
            if best is not None: break
        if best is None: best = axh.annotate(txt, xy, xytext=(6, 4), textcoords='offset points', color=col, fontsize=fs)
        boxes.append(best.get_window_extent(rd))
def meridional(T, step, fname, nper=20, skip=()):
    NAMES = [n for n in ORDER if n not in skip]
    NN = int(T*nper); TL = np.linspace(0, -T, NN + 1)
    LG = {'Sun': orbit((0.0, 0.0), ZSUN, VSUN, T=T, n=NN)}
    for p in PAST: LG[p[0]] = orbit(p[1], p[4] + ZSUN, VL[p[0]], T=T, n=NN)
    fg = plt.figure(figsize=(16, 13.5), dpi=125); fg.patch.set_facecolor(BG)
    m_ = fg.add_axes([0.07, 0.50, 0.70, 0.42]); m_.set_facecolor('#0a0f1c')
    t_ = fg.add_axes([0.07, 0.10, 0.70, 0.32]); t_.set_facecolor('#0a0f1c')
    for a_ in (m_, t_):
        for sp in a_.spines.values(): sp.set_color('#3a4a60')
        a_.tick_params(colors='#9aa', labelsize=9)
    m_.axhline(0, color='#9fb6d8', lw=.9, ls=(0, (2, 3))); m_.axvline(0, color='#3a4a60', lw=.7); t_.axhline(0, color='#3a4a60', lw=.7)
    for tk in range(-int(T), 1, step): t_.axvline(tk, color='#3a4a60', lw=.5, alpha=.6)
    PEND = []; CNT = {}
    for n in NAMES:
        px, py, pz = LG[n]; col = LCOL[n]
        CNT[n] = (int(np.sum(np.diff(np.sign(pz)) != 0)), int(np.sum(np.diff(np.sign(np.diff(py))) != 0)))
        print('%s over %d Myr: in/out %.0f to %.0f ly, height %.0f to %.0f ly' % (n, T, py.min(), py.max(), pz.min(), pz.max()))
        m_.plot(py, pz, color=col, lw=1.6 if n == 'Sun' else 1.1, alpha=.85, zorder=5)
        for tk in range(step, int(T) + 1, step):
            i = int(round(tk/T*NN)); m_.plot(py[i], pz[i], marker='o', ms=3.5, color=col, zorder=6)
            if n == 'Sun': PEND.append((m_, (py[i], pz[i]), '%d' % tk, col, 8.5))
        m_.scatter([py[-1]], [pz[-1]], s=50, facecolors='none', edgecolors=col, linewidths=1.6, zorder=7)
        m_.plot(py[0], pz[0], marker='o', ms=7, color=col, ls='none', zorder=8)
        t_.plot(TL, py, color=col, lw=2.0 if n == 'Sun' else 1.4, alpha=.95)
        t_.plot(0, py[0], marker='o', ms=7, color=col); PEND.append((t_, (0, py[0]), n, col, 10))
    m_.set_ylim(-480, 480)
    m_.set_xlabel('← toward the galactic centre        light years        away from it →', color='#ccc', fontsize=10)
    m_.set_ylabel('height above mid-plane (ly)', color='#ccc', fontsize=10)
    m_.set_title('Seen edge-on with the trip around the galaxy taken out: in-and-out against up-and-down\n(small dots every %d million years; numbers on the Sun’s path; circles = %d million years ago)' % (step, T), color='white', fontsize=11.5)
    t_.set_xlim(-T - 0.03*T, 0.18*T)
    t_.set_xlabel('millions of years from now', color='#ccc', fontsize=10); t_.set_ylabel('← toward galactic centre   (ly)   away →', color='#ccc', fontsize=10)
    t_.set_title('In and out from the galactic centre over time', color='white', fontsize=12.5)
    place_labels(fg, PEND)
    # key with the number of swings, to the right of the panels (off all lines)
    nv, nr = T*NU/(2*math.pi), T*math.sqrt(-4*OM*OORT_B*KMS/1000)/(2*math.pi)
    yk = 0.90
    fg.text(0.795, yk, 'Dot = now\nCircle = %d million years ago' % T, color='#ccc', fontsize=10, va='top'); yk -= 0.07
    for n in NAMES:
        cz, cr = CNT[n]
        fg.text(0.795, yk, n, color=LCOL[n], fontsize=12, weight='bold', va='top'); yk -= 0.028
        fg.text(0.795, yk, 'up and down: %.1f full swings\n   (crossed the mid-plane %d times)\nin and out: %.1f full swings\n   (turned around %d times)' % (nv, cz, nr, cr),
                color=LCOL[n], fontsize=10, va='top', linespacing=1.35); yk -= 0.105
    fg.text(0.795, yk, 'In this simple model every star keeps\nthe same two rhythms (%.0f and %.0f\nmillion years); only the size and timing\nof the swings differ.' % (2*math.pi/NU, T/nr),
            color='#9aa', fontsize=9.5, va='top', linespacing=1.3)
    fg.text(0.5, 0.975, 'The last %d million years, holding the motion around the galaxy constant' % T, color='white', fontsize=15, ha='center')
    fg.text(0.01, 0.005, 'Same model as pleiades_in_orion_arm.png (epicycles in the plane, Oort constants A = 15.3, B = −11.9 km/s/kpc, Bovy 2017; bounce through the mid-plane every 91 million years, McKee et al. 2015), '
            'with the drift along the orbit removed.\nZero = the Sun’s present distance from the galactic centre and the mid-plane. Over %d million years this simple model is only a sketch: spiral arms, gas clouds and the galaxy’s bar '
            'nudge real paths.\nThe Pleiades are left out: at about 120 million years old they did not exist for most of this span. The Ursa Major stars are about 400 million years old; the Hyades about 650 million.' % T, color='#9aa', fontsize=9.5)
    fg.savefig(os.path.join(out, fname), facecolor=BG, bbox_inches='tight'); print('saved', fname)
meridional(360, 40, 'pleiades_360_million_years_in_out_up_down.png', skip=('Pleiades',))

# ---- seventh figure: paths around the galaxy over the next billion years (rosettes) for the Sun, Pleiades, Hyades and Ursa Major core ----
TF = 1000.0; NF = 20000
tf = np.linspace(0, TF, NF + 1)
R0LY = R0*KPC; OM0 = OM
kap = math.sqrt(-4*OM*OORT_B*KMS/1000)
ROS = {}
for n in ORDER:
    if n == 'Sun': ux, uy, uz = orbit((0.0, 0.0), ZSUN, VSUN, T=-TF, n=NF)
    else:
        p = [q for q in PAST if q[0] == n][0]; ux, uy, uz = orbit(p[1], p[4] + ZSUN, VL[n], T=-TF, n=NF)
    th = math.pi/2 - OM0*tf - ux/R0LY; rr = R0LY + uy
    apo = [i for i in range(1, NF) if rr[i] > rr[i - 1] and rr[i] >= rr[i + 1]]
    steps = [math.degrees(th[apo[k]] - th[apo[k + 1]]) for k in range(len(apo) - 1)]
    ROS[n] = (uy, th, apo, uy.mean())
    print('%s: average distance %+.0f ly from the circle; farthest-out points at %s Myr; each %.0f deg further around (%.0f deg short of a full turn)'
          % (n, uy.mean(), [int(round(tf[a], -1)) for a in apo], np.mean(steps), 360 - np.mean(steps)))
fig7 = plt.figure(figsize=(19, 9.4), dpi=125); fig7.patch.set_facecolor(BG)
for k, (ex, tt) in enumerate([(1, 'True scale'), (3, 'Same paths with each in-and-out swing drawn 3 times larger')]):
    rx = fig7.add_axes([0.02 + 0.40*k, 0.08, 0.38, 0.76]); rx.set_facecolor('#0a0f1c')
    for sp in rx.spines.values(): sp.set_color('#3a4a60')
    rx.set_xticks([]); rx.set_yticks([])
    cc = np.linspace(0, 2*math.pi, 400); rx.plot(R0LY*np.cos(cc), R0LY*np.sin(cc), color='#ffd36e', lw=.9, ls=(0, (5, 4)), alpha=.6)
    rx.plot(0, 0, marker='*', ms=14, color='#ffb347')
    pend = [(rx, (0, 0), 'galactic centre', '#ffb347', 10)]; lim = 0
    for n in ORDER:
        uy, th, apo, ug = ROS[n]; col = LCOL[n]
        Rp = R0LY + ug + ex*(uy - ug); X, Y = Rp*np.cos(th), Rp*np.sin(th); lim = max(lim, Rp.max())
        rx.plot(X, Y, color=col, lw=1.5 if n == 'Sun' else 0.9, alpha=.85, zorder=5 if n == 'Sun' else 4)
        for a in apo: rx.plot(X[a], Y[a], marker='o', ms=6, mfc=BG, mec=col, mew=1.4, zorder=7)
        if n == 'Sun':
            for a in apo: pend.append((rx, (X[a], Y[a]), '%d' % round(tf[a], -1), col, 9))
        rx.plot(X[0], Y[0], marker='o', ms=8, color=col, zorder=8)
    lim *= 1.07; rx.set_xlim(-lim, lim); rx.set_ylim(-lim, lim); rx.set_aspect('equal')
    rx.set_title(tt, color='white', fontsize=12.5)
    place_labels(fig7, pend)
yk = 0.84
fig7.text(0.815, yk, 'Filled dot = now (all four start\nwithin 450 ly of each other, at the top)\nRing = a farthest-out point', color='#ccc', fontsize=10, va='top', linespacing=1.3); yk -= 0.10
for n in ORDER:
    uy, th, apo, ug = ROS[n]
    steps = [math.degrees(th[apo[k]] - th[apo[k + 1]]) for k in range(len(apo) - 1)]
    fig7.text(0.815, yk, n, color=LCOL[n], fontsize=12, weight='bold', va='top'); yk -= 0.03
    fig7.text(0.815, yk, 'on average %s ly %s the circle\nfarthest out every 171 million years,\neach time about %.0f° short of a full turn' %
              (format(int(round(abs(ug), -1)), ','), 'outside' if ug > 0 else 'inside', 360 - np.mean(steps)), color=LCOL[n], fontsize=10, va='top', linespacing=1.3); yk -= 0.115
fig7.text(0.815, yk, 'The Sun’s ring labels are millions\nof years from now.', color='#9aa', fontsize=9.5, va='top', linespacing=1.3)
fig7.text(0.42, 0.955, 'Paths around the galaxy over the next billion years, seen from above', color='white', fontsize=15, ha='center')
fig7.text(0.42, 0.915, 'None of the paths closes on itself: each farthest-out point lands short of the last, so every path traces a rosette.', color='#e8e8e8', fontsize=11, ha='center')
fig7.text(0.01, 0.01, 'Dashed circle = a circular orbit at the Sun’s distance (%s ly from the centre; one turn every %.0f million years; the galaxy turns clockwise in this view). '
          'Same simple model as the other pictures (epicycles: Oort constants A = 15.3, B = −11.9 km/s/kpc, Bovy 2017), run forward.\nThe stars that drift outward fall behind and those that drift inward move ahead, so the four spread around the galaxy. '
          'Over a billion years spiral arms, the bar and gas clouds will also nudge the real paths, and the clusters will have dispersed: this shows the pattern, not a forecast.'
          % (format(int(round(R0LY, -2)), ','), 2*math.pi/OM0), color='#9aa', fontsize=9.5)
fig7.savefig(os.path.join(out, 'rosettes_1_billion_years.png'), facecolor=BG, bbox_inches='tight'); print('saved rosette')
