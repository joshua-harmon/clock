# Sirius A's path over 500 years including the pull of Sirius B: the barycentre of the pair moves in a (nearly)
# straight line with the Hipparcos proper motion, and Sirius A loops around it once every 50.13 years.
# Orbit of Sirius B: Bond et al. 2017 (ApJ 840, 70); masses A 2.063, B 1.018 solar masses.
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sirius_orbit import a_off, rel, P, F, A as AREL

here = os.path.dirname(os.path.abspath(__file__))
K = 4.740470446; KMS_LY = 1.0227e-6*3.26156
RA0, DE0, PLX, PMRA, PMDE, RV = 101.2872, -16.7161, 379.21, -546.01, -1223.07, -5.5; EP0 = 1991.25
def unit(ra, de):
    a, d = math.radians(ra), math.radians(de); return np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
r0 = unit(RA0, DE0); a, d = math.radians(RA0), math.radians(DE0)
p = np.array([-math.sin(a), math.cos(a), 0]); q = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
D0 = 3261.56/PLX; vel = ((K*PMRA/PLX)*p + (K*PMDE/PLX)*q + RV*r0)*KMS_LY
def bary(y):            # barycentre, arcsec east/north of a fixed reference
    v = D0*r0 + (y - EP0)*vel; u = v/np.linalg.norm(v); return 3600*math.degrees(1)*np.array([np.dot(u, p), np.dot(u, q)])
def siriusA(y): return bary(y) + a_off(y)
REF = siriusA(2026)
def A_(y): return (siriusA(y) - REF)/60          # arcmin from Sirius A's 2026 place
def B_(y): return (bary(y) - REF)/60
yrs = np.arange(1526, 2026.001, 0.1)
sa = np.array([A_(y) for y in yrs]); sb = np.array([B_(y) for y in yrs])
dirn = np.array([PMRA, PMDE]); dirn = dirn/np.linalg.norm(dirn); perp = np.array([-dirn[1], dirn[0]])
offs = np.array([a_off(y) for y in yrs])        # arcsec
along, across = offs @ dirn, offs @ perp
print('max wobble %.2f" ; across-track range %.2f" ; along-track range %.2f"' % (np.max(np.hypot(*offs.T)), np.ptp(across), np.ptp(along)))

BG = '#05070d'; HIST = [(1718, 'Halley'), (1844, 'Bessel'), (1862, 'Clark sees B')]
def style(ax, title):
    ax.set_facecolor(BG); ax.tick_params(colors='#9aa', labelsize=9); ax.grid(color='#1c2636', lw=.6)
    for s_ in ax.spines.values(): s_.set_color('#3a4a60')
    ax.set_title(title, color='white', fontsize=12)

fig = plt.figure(figsize=(21, 13), dpi=120); fig.patch.set_facecolor(BG)
ax1 = fig.add_axes([0.04, 0.09, 0.40, 0.80]); ax2 = fig.add_axes([0.50, 0.53, 0.22, 0.36]); ax3 = fig.add_axes([0.76, 0.53, 0.22, 0.36]); ax4 = fig.add_axes([0.50, 0.09, 0.48, 0.33])

# 1. the whole 500 years on the sky
style(ax1, 'Sirius A, 1526 to 2026, with the pull of Sirius B')
ax1.plot(sb[:, 0], sb[:, 1], color='#888', lw=1, ls=(0, (4, 3)), label='centre of mass of Sirius A and B (straight)')
ax1.plot(sa[:, 0], sa[:, 1], color='#7fb2ff', lw=1.5, label='Sirius A (loops around it every 50.13 years)')
for y in range(1550, 2026, 50):
    x_, y_ = A_(y); ax1.scatter([x_], [y_], s=36, color='#e6ecff', edgecolors='#7fb2ff', zorder=3)
    ax1.annotate(str(y), (x_, y_), xytext=(-10, 0), textcoords='offset points', color='#9fc4ff', fontsize=10, ha='right', va='center')
for y, nm in HIST:
    x_, y_ = A_(y); ax1.scatter([x_], [y_], s=55, color='#ffd27a', zorder=4)
    ax1.annotate('%d %s' % (y, nm), (x_, y_), xytext=(10, 0), textcoords='offset points', color='#ffd27a', fontsize=10, va='center')
ax1.scatter([0], [0], s=240, color='white', zorder=5); ax1.annotate('Sirius A in 2026', (0, 0), xytext=(14, -2), textcoords='offset points', color='white', fontsize=12, weight='bold', va='center')
ax1.add_patch(plt.Rectangle((0.5, 3.0), 2.6, -2.25, fill=False, ec='#ffd27a', lw=.8, ls=':'))   # placeholder, replaced below
ax1.patches[-1].remove()
x0, y0 = A_(1830); x1, y1 = A_(1880)
ax1.add_patch(plt.Rectangle((min(x0, x1) - .15, min(y0, y1) - .15), abs(x1 - x0) + .3, abs(y1 - y0) + .3, fill=False, ec='#c7a6ff', lw=1))
ax1.text(min(x0, x1) - .25, max(y0, y1) + .1, 'enlarged at\ntop middle', color='#c7a6ff', fontsize=9, ha='left')
ax1.set_xlim(5.6, -5.0); ax1.set_ylim(-0.9, 11.0); ax1.set_aspect('equal')
ax1.set_xlabel('arcminutes east of Sirius A’s 2026 position (east to the left)', color='#ccc'); ax1.set_ylabel('arcminutes north', color='#ccc')
ax1.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=9.5)
ax1.text(-4.8, 10.6, 'The wobble is at most 3.7″ (0.06′):\nsmall at this scale, see the panels at right', color='#9aa', fontsize=9.5, ha='right', va='top')

# 2. close-up 1830-1880 (Bessel's years), same scale both ways, in arcseconds about the 1855 barycentre
style(ax2, 'Close-up, 1830 to 1880 (one orbit), in arcseconds')
yy = np.arange(1830, 1880.001, 0.1); c = bary(1855)
pa = np.array([siriusA(y) - c for y in yy]); pb = np.array([bary(y) - c for y in yy])
ax2.plot(pb[:, 0], pb[:, 1], color='#888', lw=1, ls=(0, (4, 3))); ax2.plot(pa[:, 0], pa[:, 1], color='#7fb2ff', lw=1.6)
for y in range(1830, 1881, 10):
    v = siriusA(y) - c; ax2.scatter([v[0]], [v[1]], s=26, color='#e6ecff', zorder=3)
    ax2.annotate(str(y), v, xytext=(-8, 0), textcoords='offset points', color='#9fc4ff', fontsize=9, ha='right', va='center')
for y in (1844, 1862):
    v = siriusA(y) - c; ax2.scatter([v[0]], [v[1]], s=45, color='#ffd27a', zorder=4)
    ax2.annotate(str(y), v, xytext=(8, 0), textcoords='offset points', color='#ffd27a', fontsize=9, va='center')
ax2.set_aspect('equal'); ax2.invert_xaxis(); ax2.set_xlabel('arcseconds east', color='#ccc'); ax2.set_ylabel('arcseconds north', color='#ccc')

# 3. Sirius A around the centre of mass, and Sirius B around it, arcseconds
style(ax3, 'The two orbits about the centre of mass')
ty = np.arange(0, P, 0.05) + 1994.5715
oa = np.array([a_off(t) for t in ty]); ob = np.array([a_off(t) + rel(t) for t in ty])
ax3.plot(ob[:, 0], ob[:, 1], color='#bbbbbb', lw=1.1, label='Sirius B'); ax3.plot(oa[:, 0], oa[:, 1], color='#7fb2ff', lw=1.8, label='Sirius A')
ax3.scatter([0], [0], marker='+', s=80, color='#ffe24a', zorder=4); ax3.annotate('centre of mass', (0, 0), xytext=(-4, 10), textcoords='offset points', color='#ffe24a', fontsize=8.5)
for y in (1844, 1862, 1994.6, 2026):
    va, vb = a_off(y), a_off(y) + rel(y)
    ax3.plot([va[0], vb[0]], [va[1], vb[1]], color='#ffd27a' if y < 1900 else '#556', lw=.6, ls=':')
    lab = {1994.6: '', 1844: '1844 and 1994 (closest)'}.get(y, str(int(y)))
    ax3.scatter([va[0]], [va[1]], s=30, color='#ffd27a' if y < 1900 else 'white', zorder=4); ax3.scatter([vb[0]], [vb[1]], s=14, color='#ccc', zorder=4)
    ax3.annotate(lab, vb, xytext=(5, -12) if y == 1844 else (5, 3), textcoords='offset points', color='#ddd', fontsize=8.5)
ax3.set_aspect('equal'); ax3.invert_xaxis(); ax3.set_xlabel('arcseconds east', color='#ccc'); ax3.set_ylabel('arcseconds north', color='#ccc')
ax3.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=8.5)

# 4. the wobble over 500 years, split along and across the direction of travel
style(ax4, 'Sirius A’s wobble over 500 years: its offset from the steady (centre-of-mass) motion')
ax4.plot(yrs, along, color='#7fb2ff', lw=1.6, label='ahead (+) or behind (−) along the direction of travel')
ax4.plot(yrs, across, color='#ff9f6b', lw=1.4, label='to the side of the track')
ax4.axhline(0, color='#556', lw=.8)
for (y, nm), ty_ in zip(HIST, (3.9, 3.9, 4.35)): ax4.axvline(y, color='#ffd27a', lw=.7, ls=':'); ax4.text(y + 2, ty_, '%d %s' % (y, nm), color='#ffd27a', fontsize=9)
ax4.set_xlim(1526, 2026); ax4.set_ylim(-4.3, 4.9); ax4.set_xlabel('year', color='#ccc'); ax4.set_ylabel('arcseconds', color='#ccc')
ax4.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=9, ncol=2)

fig.suptitle('Sirius over 500 years, including the pull of its white-dwarf companion Sirius B', color='white', fontsize=15.5)
fig.text(.04, .015, 'The centre of mass moves steadily at the Hipparcos (2007) proper motion, 1.339″ a year toward the south-southwest (Hipparcos fitted this binary with its orbit, so its proper motion is the steady one). '
         'Sirius B’s orbit: period 50.13 years, 7.50″ across the semi-major axis, eccentricity 0.59 (Bond et al. 2017); Sirius A (2.06 solar masses) swings 33% as far as Sirius B (1.02) on the opposite side. '
         'Fixed J2000 coordinates, no precession. In 1844 Bessel found from a century of positions that Sirius sped up and slowed down along its track, the along-track wobble in the bottom panel.',
         color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, 'sirius_motion_500_years_with_sirius_b.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print('done')
