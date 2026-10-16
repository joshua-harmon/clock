# Sky chart (500 years with the companion's pull) and north-galactic-pole chart for one binary: python make_binary_charts.py Procyon
import math, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from binary_stars import Binary, top

here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
s = Binary(sys.argv[1] if len(sys.argv) > 1 else 'Procyon'); N = s.name; c0, c1 = s.close
def style(ax, title):
    ax.set_facecolor(BG); ax.tick_params(colors='#9aa', labelsize=9); ax.grid(color='#1c2636', lw=.6)
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.set_title(title, color='white', fontsize=12)

# ---------------- chart 1: on the sky ----------------
REF = s.sky_A(2026)
def A_(y): return (s.sky_A(y) - REF)/60
def B_(y): return (s.sky_bary(y) - REF)/60
yrs = np.arange(1526, 2026.001, 0.1)
sa = np.array([A_(y) for y in yrs]); sb = np.array([B_(y) for y in yrs])
dirn = np.array([s.pmra, s.pmde]); dirn /= np.linalg.norm(dirn); perp = np.array([-dirn[1], dirn[0]])
offs = np.array([s.a_off(y) for y in yrs]); along, across = offs @ dirn, offs @ perp
tot = math.hypot(*A_(1526)); wob = np.max(np.hypot(*offs.T))

fig = plt.figure(figsize=(21, 13), dpi=120); fig.patch.set_facecolor(BG)
ax1 = fig.add_axes([0.04, 0.09, 0.40, 0.80]); ax2 = fig.add_axes([0.50, 0.53, 0.22, 0.36]); ax3 = fig.add_axes([0.76, 0.53, 0.22, 0.36]); ax4 = fig.add_axes([0.50, 0.09, 0.48, 0.33])
style(ax1, '%s A, 1526 to 2026, with the pull of %s' % (N, s.B))
ax1.plot(sb[:, 0], sb[:, 1], color='#888', lw=1, ls=(0, (4, 3)), label='centre of mass of %s A and B (straight)' % N)
ax1.plot(sa[:, 0], sa[:, 1], color=s.colour, lw=1.5, label='%s A (loops around it every %.2f years)' % (N, s.P))
for y in range(1550, 2026, 50):
    x_, y_ = A_(y); ax1.scatter([x_], [y_], s=36, color='#e6ecff', edgecolors=s.colour, zorder=3)
    ax1.annotate(str(y), (x_, y_), xytext=(-10, 0), textcoords='offset points', color='#cfd8e6', fontsize=10, ha='right', va='center')
for k, (y, nm) in enumerate(s.hist):
    x_, y_ = A_(y); ax1.scatter([x_], [y_], s=55, color='#ffb36b', zorder=4)
    ax1.annotate('%d %s' % (y, nm), (x_, y_), xytext=[(12, 0), (12, 12), (12, -3), (12, -8)][k], textcoords='offset points', color='#ffb36b', fontsize=10, va='center')
ax1.scatter([0], [0], s=240, color='white', zorder=5); ax1.annotate('%s A in 2026' % N, (0, 0), xytext=(14, -2), textcoords='offset points', color='white', fontsize=12, weight='bold', va='center')
x0, y0 = A_(c0); x1, y1 = A_(c1)
ax1.add_patch(plt.Rectangle((min(x0, x1) - .15, min(y0, y1) - .15), abs(x1 - x0) + .3, abs(y1 - y0) + .3, fill=False, ec='#c7a6ff', lw=1))
ax1.plot([], [], color='#c7a6ff', lw=0, marker='s', mfc='none', ms=10, label='box: enlarged at top middle')
xs, ys = sa[:, 0], sa[:, 1]
cx_, cy_ = (xs.max() + xs.min())/2, (ys.max() + ys.min())/2; half = max(np.ptp(xs), np.ptp(ys))/2 + 1.0
ax1.set_xlim(cx_ + half, cx_ - half); ax1.set_ylim(cy_ - half, cy_ + half); ax1.set_aspect('equal')
ax1.set_xlabel('arcminutes east of %s A’s 2026 position (east to the left)' % N, color='#ccc'); ax1.set_ylabel('arcminutes north', color='#ccc')
ax1.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=9.5)
ax1.text(cx_ - half + .2, cy_ + half - .2, 'The wobble is at most %.1f″ (%.2f′):\nsmall at this scale, see the panels at right' % (wob, wob/60), color='#9aa', fontsize=9.5, ha='right', va='top')

style(ax2, 'Close-up, %d to %d (one orbit), in arcseconds' % (c0, c1))
yy = np.arange(c0, c1 + .001, 0.1); cm = s.sky_bary((c0 + c1)/2)
pa = np.array([s.sky_A(y) - cm for y in yy]); pb = np.array([s.sky_bary(y) - cm for y in yy])
ax2.plot(pb[:, 0], pb[:, 1], color='#888', lw=1, ls=(0, (4, 3))); ax2.plot(pa[:, 0], pa[:, 1], color=s.colour, lw=1.6)
for y in range(c0 - c0 % 10 + 10, c1 + 1, 10):
    v = s.sky_A(y) - cm; ax2.scatter([v[0]], [v[1]], s=26, color='#e6ecff', zorder=3)
    ax2.annotate(str(y), v, xytext=(-8, 0), textcoords='offset points', color='#cfd8e6', fontsize=9, ha='right', va='center')
for y, nm in s.hist:
    if c0 <= y <= c1:
        v = s.sky_A(y) - cm; ax2.scatter([v[0]], [v[1]], s=45, color='#ffb36b', zorder=4)
        ax2.annotate(str(y), v, xytext=(8, 0), textcoords='offset points', color='#ffb36b', fontsize=9, va='center')
ax2.set_aspect('equal'); ax2.invert_xaxis(); ax2.set_xlabel('arcseconds east', color='#ccc'); ax2.set_ylabel('arcseconds north', color='#ccc')

style(ax3, 'The two orbits about the centre of mass (sky)')
ty = np.arange(0, s.P, 0.05) + s.T0
oa = np.array([s.a_off(t) for t in ty]); ob = np.array([s.a_off(t) + s.rel(t)[:2] for t in ty])
ax3.plot(ob[:, 0], ob[:, 1], color='#bbbbbb', lw=1.1, label=s.B); ax3.plot(oa[:, 0], oa[:, 1], color=s.colour, lw=1.8, label='%s A' % N)
ax3.scatter([0], [0], marker='+', s=80, color='#ffe24a', zorder=4)
for y in [h[0] for h in s.hist[1:]] + [2026]:
    va, vb = s.a_off(y), s.a_off(y) + s.rel(y)[:2]; col = 'white' if y == 2026 else '#ffb36b'
    ax3.plot([va[0], vb[0]], [va[1], vb[1]], color=col, lw=.6, ls=':')
    ax3.scatter([va[0]], [va[1]], s=30, color=col, zorder=4); ax3.scatter([vb[0]], [vb[1]], s=14, color=col, zorder=4)
    ax3.annotate(str(y), vb, xytext=(5, -11) if y == 1862 else (5, 3), textcoords='offset points', color='#ddd', fontsize=8.5)
ax3.set_aspect('equal'); ax3.invert_xaxis(); ax3.margins(0.15); ax3.set_xlabel('arcseconds east', color='#ccc'); ax3.set_ylabel('arcseconds north', color='#ccc')
ax3.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=8.5)

style(ax4, '%s A’s wobble over 500 years: its offset from the steady (centre-of-mass) motion' % N)
ax4.plot(yrs, along, color=s.colour, lw=1.6, label='ahead (+) or behind (−) along the direction of travel')
ax4.plot(yrs, across, color='#ff9f6b', lw=1.4, label='to the side of the track')
ax4.axhline(0, color='#556', lw=.8); lim = max(np.abs(along).max(), np.abs(across).max())
for k, (y, nm) in enumerate(s.hist): ax4.axvline(y, color='#ffb36b', lw=.7, ls=':'); ax4.text(y + 2, lim*(1.18 + .12*(k % 2)), str(y), color='#ffb36b', fontsize=9)
ax4.set_xlim(1526, 2026); ax4.set_ylim(-lim*1.45, lim*1.45); ax4.set_xlabel('year', color='#ccc'); ax4.set_ylabel('arcseconds', color='#ccc')
ax4.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=9, ncol=2)
fig.suptitle('%s over 500 years, including the pull of its white-dwarf companion %s: %.1f arcminutes in all' % (N, s.B, tot), color='white', fontsize=15.5)
fig.text(.04, .015, 'Centre of mass: Hipparcos (2007) proper motion, %.3f″ a year toward position angle %.0f°, radial velocity %.1f km/s; taken as the steady motion (Hipparcos fitted these binaries with their known orbits). '
         '%s’s orbit: period %.2f years, semi-major axis %.2f″, eccentricity %.2f (Bond et al.); %s A (%.2f suns) swings %.0f%% as far as %s (%.2f suns) on the opposite side. Fixed J2000 coordinates, no precession.'
         % (math.hypot(s.pmra, s.pmde)/1000, math.degrees(math.atan2(s.pmra, s.pmde)) % 360, s.rv, s.B, s.P, s.a, s.e, N, s.mA, 100*s.F/(1 - s.F), s.B, s.mB),
         color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, '%s_motion_500_years_with_%s_b.png' % (N.lower(), N.lower())), facecolor=BG, bbox_inches='tight'); plt.close(fig)

# ---------------- chart 2: from the north galactic pole ----------------
sp, vin, vz = s.speeds(); D = s.D0
fig = plt.figure(figsize=(21, 7.8 if N == 'Procyon' else 11.5), dpi=120); fig.patch.set_facecolor(BG)
ax1 = fig.add_axes([0.03, 0.10, 0.22, 0.78]); axo = fig.add_axes([0.31, 0.10, 0.40, 0.78]); ax2 = fig.add_axes([0.77, 0.10, 0.21, 0.78])
C = top(s.bary3((c0 + c1)/2)); yy = np.arange(c0, c1 + .001, 0.05)
cb = np.array([top(s.bary3(t)) - C for t in yy]); ca = np.array([top(s.A3(t)) - C for t in yy]); cB = np.array([top(s.B3(t)) - C for t in yy])
style(ax1, 'Path through space, %d to %d' % (c0, c1))
ax1.plot(cb[:, 0], cb[:, 1], color='#888', lw=1, ls=(0, (4, 3))); ax1.plot(cB[:, 0], cB[:, 1], color='#c8c8c8', lw=1.1); ax1.plot(ca[:, 0], ca[:, 1], color=s.colour, lw=2)
for t in range(c0 - c0 % 10 + 10, c1 + 1, 10):
    va = top(s.A3(t)) - C; ax1.scatter([va[0]], [va[1]], s=50, color='#e6ecff', edgecolors=s.colour, zorder=4)
    ax1.annotate(str(t), va, xytext=(-9, 0), textcoords='offset points', color='#cfd8e6', fontsize=10, ha='right', va='center')
for y, nm in s.hist:
    if c0 <= y <= c1:
        va = top(s.A3(y)) - C; ax1.scatter([va[0]], [va[1]], s=70, color='#ffb36b', zorder=5)
        ax1.annotate(str(y), va, xytext=(9, 0), textcoords='offset points', color='#ffb36b', fontsize=10, va='center')
ax1.set_aspect('equal'); ax1.margins(0.12)
ax1.set_xlabel('AU toward l = 90° (Cygnus) →', color='#ccc'); ax1.set_ylabel('AU  (galactic centre is down)', color='#ccc')

style(axo, 'Riding along with the pair: both orbits seen from the north galactic pole')
to = np.arange(0, s.P, 0.02) + s.T0
oa = np.array([top(s.A3(t) - s.bary3(t)) for t in to]); ob = np.array([top(s.B3(t) - s.bary3(t)) for t in to])
axo.plot(ob[:, 0], ob[:, 1], color='#c8c8c8', lw=1.2, label='%s (white dwarf, %.2f suns)' % (s.B, s.mB))
axo.plot(oa[:, 0], oa[:, 1], color=s.colour, lw=2.2, label='%s A (%.2f suns)' % (N, s.mA))
axo.scatter([0], [0], marker='+', s=120, color='#ffe24a', zorder=5, label='centre of mass')
peri = s.T0 + s.P*math.ceil((1830 - s.T0)/s.P); apo = peri + s.P/2
R = 1.15*max(np.abs(ob).max(), 1)
marks = [(peri, '%d and %d: closest, %.1f AU apart' % (round(peri), round(peri + s.P*math.floor((2026 - peri)/s.P)), s.a*s.au*(1 - s.e)), '#bbb'),
         (apo, '%d and %d: farthest, %.1f AU apart' % (round(apo), round(apo + s.P*math.floor((2026 - apo)/s.P)), s.a*s.au*(1 + s.e)), '#bbb')]
marks += [(y, '%d: %s' % (y, nm), '#ffb36b') for y, nm in s.hist[1:]] + [(2026, '2026', 'white')]
marks.sort(key=lambda m: -top(s.B3(m[0]) - s.bary3(m[0]))[1]); ys_lab = np.linspace(0.85*R, -0.85*R, len(marks))
for (t, lab, col), ty in zip(marks, ys_lab):
    va, vb = top(s.A3(t) - s.bary3(t)), top(s.B3(t) - s.bary3(t))
    axo.plot([va[0], vb[0]], [va[1], vb[1]], color=col, lw=.8, ls=':')
    axo.scatter([va[0]], [va[1]], s=55, color=col, zorder=4); axo.scatter([vb[0]], [vb[1]], s=25, color=col, zorder=4)
    axo.annotate(lab, vb, xytext=(1.25*R, ty), textcoords='data', color=col, fontsize=10, va='center', arrowprops=dict(arrowstyle='-', color=col, lw=.6, alpha=.5, shrinkB=4))
for t in np.arange(s.T0, s.T0 + s.P, 5):
    vb = top(s.B3(t) - s.bary3(t)); axo.scatter([vb[0]], [vb[1]], s=8, color='#888', zorder=3)
axo.set_xlim(-1.3*R, 2.9*R); axo.set_ylim(-1.15*R, 1.15*R); axo.set_aspect('equal')
axo.set_xlabel('AU toward galactic longitude 90° (Cygnus)  →', color='#ccc'); axo.set_ylabel('AU  (galactic centre is down)', color='#ccc')
axo.legend(loc='upper left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=10)
axo.text(-1.25*R, -1.08*R, 'Dotted lines join %s A (big dot) to %s (small dot) at that date.\nSmall grey dots on B’s orbit: every 5 years.' % (N, s.B), color='#9aa', fontsize=9)
sun = top(-s.bary3(1900)); sun /= np.linalg.norm(sun); mv = top(np.array(__import__('binary_stars').T @ s.vel)); mv /= np.linalg.norm(mv)
cx, cy, L = -0.85*R, -0.45*R, 0.22*R
for vec, col, txt, ha in ((sun, '#ffe24a', 'toward the Sun', 'left'), (np.array([0, -1.0]), '#ff9f6b', 'galactic centre', 'center'), (mv, s.colour, 'travel', 'right')):
    axo.annotate('', xy=(cx + L*vec[0], cy + L*vec[1]), xytext=(cx, cy), arrowprops=dict(arrowstyle='-|>', color=col, lw=1.8))
    axo.text(cx + 1.25*L*vec[0], cy + 1.25*L*vec[1], txt, color=col, fontsize=9, ha=ha, va='center')

C2 = top(s.bary3(2026)); y5 = np.arange(1526, 2026.001, 0.2)
wa = np.array([top(s.A3(t)) - C2 for t in y5]); wb = np.array([top(s.bary3(t)) - C2 for t in y5])
style(ax2, 'The whole 500 years, same view')
ax2.plot(wb[:, 0], wb[:, 1], color='#888', lw=1, ls=(0, (4, 3))); ax2.plot(wa[:, 0], wa[:, 1], color=s.colour, lw=1.4)
for t in range(1550, 2026, 50):
    v = top(s.A3(t)) - C2; ax2.scatter([v[0]], [v[1]], s=26, color='#e6ecff', zorder=3)
    ax2.annotate(str(t), v, xytext=(8, 0), textcoords='offset points', color='#cfd8e6', fontsize=9, va='center')
v = top(s.A3(2026)) - C2; ax2.scatter([v[0]], [v[1]], s=160, color='white', zorder=4); ax2.annotate('2026', v, xytext=(9, -8), textcoords='offset points', color='white', fontsize=10.5, weight='bold')
ax2.set_aspect('equal'); ax2.margins(0.15); ax2.set_xlabel('AU toward l = 90°  →', color='#ccc'); ax2.set_ylabel('AU (galactic centre down)', color='#ccc')
fig.suptitle('%s A and B seen from the north galactic pole: their real path through space relative to the Sun' % N, color='white', fontsize=15)
fig.text(.04, .015, 'Top view: looking down on the Milky Way’s disc from the north, Cygnus (l = 90°) to the right, the galactic centre down. %s is %.1f ly away and moves at %.1f km/s relative to the Sun '
         '(%.1f km/s along the plane, %.1f km/s toward the %s), so the pair advances about %d AU across the plane in each %.1f-year orbit. 1 AU = the Earth–Sun distance. '
         'The front/back sense of the orbit’s tilt is taken from the usual radial-velocity convention; if reversed, the middle panel would be mirrored.'
         % (N, D, sp, vin, abs(vz), 'south' if vz < 0 else 'north', round(vin*s.P*1.0227e-6*3.26156*63241.077, -1), s.P), color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, '%s_motion_from_north_galactic_pole.png' % N.lower()), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print(N, 'total %.2f arcmin, wobble %.2f"' % (tot, wob), 'along %.2f across %.2f' % (np.ptp(along), np.ptp(across)))
