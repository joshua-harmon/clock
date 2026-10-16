# Sirius A and B in space, seen from the north galactic pole (top view: l = 90° to the right, galactic centre down).
# Barycentre: Hipparcos (2007) position, parallax, proper motion and RV -5.5 km/s, straight-line motion relative to the Sun.
# Orbit: Bond et al. 2017; the node is taken as the one where Sirius B recedes from us (the usual convention when
# radial velocities fix it).
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sirius_orbit import P, A as AREL, E, I, NODE, OM, T0, F

here = os.path.dirname(os.path.abspath(__file__))
T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]])
K = 4.740470446; KMS_LY = 1.0227e-6*3.26156; AU_LY = 63241.077
RA0, DE0, PLX, PMRA, PMDE, RV = 101.2872, -16.7161, 379.21, -546.01, -1223.07, -5.5; EP0 = 1991.25
a, d = math.radians(RA0), math.radians(DE0)
r0 = np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
p = np.array([-math.sin(a), math.cos(a), 0]); q = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
D0 = 3261.56/PLX; vel = ((K*PMRA/PLX)*p + (K*PMDE/PLX)*q + RV*r0)*KMS_LY
AU_PER_ARCSEC = 1/(PLX/1000)
def rel3(t):   # B relative to A, AU, equatorial frame
    M = 2*math.pi*(t - T0)/P; Ea = M
    for _ in range(60): Ea -= (Ea - E*math.sin(Ea) - M)/(1 - E*math.cos(Ea))
    nu = 2*math.atan2(math.sqrt(1+E)*math.sin(Ea/2), math.sqrt(1-E)*math.cos(Ea/2))
    r = AREL*(1 - E*math.cos(Ea))*AU_PER_ARCSEC; u = math.radians(OM) + nu; O, i = math.radians(NODE), math.radians(I)
    north = r*(math.cos(O)*math.cos(u) - math.sin(O)*math.sin(u)*math.cos(i))
    east = r*(math.sin(O)*math.cos(u) + math.cos(O)*math.sin(u)*math.cos(i))
    away = r*math.sin(u)*math.sin(i)
    return east*p + north*q + away*r0
def bary(t): return T @ ((D0*r0 + (t - EP0)*vel)*AU_LY)           # galactic, AU from the Sun
def posA(t): return bary(t) - F*(T @ rel3(t))
def posB(t): return bary(t) + (1 - F)*(T @ rel3(t))
def top(v): return np.array([v[1], -v[0]])                           # right = toward l = 90°, down = toward the galactic centre

g = T @ vel; vt = math.hypot(g[0], g[1])/KMS_LY; print('speed %.1f km/s, in-plane %.1f km/s, vertical %.1f km/s' % (np.linalg.norm(vel)/KMS_LY, vt, g[2]/KMS_LY))
print('in-plane AU per 50 yr: %.0f' % (vt*KMS_LY*50*AU_LY))
BG = '#05070d'
def style(ax, title):
    ax.set_facecolor(BG); ax.tick_params(colors='#9aa', labelsize=9); ax.grid(color='#1c2636', lw=.6)
    for s_ in ax.spines.values(): s_.set_color('#3a4a60')
    ax.set_title(title, color='white', fontsize=12.5)

fig = plt.figure(figsize=(21, 11.5), dpi=120); fig.patch.set_facecolor(BG)
ax1 = fig.add_axes([0.03, 0.10, 0.22, 0.78]); axo = fig.add_axes([0.31, 0.10, 0.40, 0.78]); ax2 = fig.add_axes([0.77, 0.10, 0.21, 0.78])

# 1. close-up 1830-1880 path, AU about the 1855 centre of mass
C = top(bary(1855)); yy = np.arange(1830, 1880.001, 0.05)
cb = np.array([top(bary(t)) - C for t in yy]); ca = np.array([top(posA(t)) - C for t in yy]); cB = np.array([top(posB(t)) - C for t in yy])
style(ax1, 'Path through space, 1830 to 1880')
ax1.plot(cb[:, 0], cb[:, 1], color='#888', lw=1, ls=(0, (4, 3)))
ax1.plot(cB[:, 0], cB[:, 1], color='#c8c8c8', lw=1.1); ax1.plot(ca[:, 0], ca[:, 1], color='#7fb2ff', lw=2)
for t in range(1830, 1881, 10):
    va = top(posA(t)) - C
    ax1.scatter([va[0]], [va[1]], s=50, color='#e6ecff', edgecolors='#7fb2ff', zorder=4)
    ax1.annotate(str(t), va, xytext=(-9, 0), textcoords='offset points', color='#9fc4ff', fontsize=10, ha='right', va='center')
for t, txt in ((1844, '1844'), (1862, '1862')):
    va = top(posA(t)) - C; ax1.scatter([va[0]], [va[1]], s=70, color='#ffd27a', zorder=5)
    ax1.annotate(txt, va, xytext=(9, 0), textcoords='offset points', color='#ffd27a', fontsize=10, va='center')
ax1.set_aspect('equal'); ax1.margins(0.12)
ax1.set_xlabel('AU toward l = 90° (Cygnus) →', color='#ccc'); ax1.set_ylabel('AU  (galactic centre is down)', color='#ccc')

# 2. the orbits with the steady motion taken out (riding along with the centre of mass)
style(axo, 'Riding along with the pair: both orbits seen from the north galactic pole')
to = np.arange(0, P, 0.02) + T0
oa = np.array([top(posA(t) - bary(t)) for t in to]); ob = np.array([top(posB(t) - bary(t)) for t in to])
axo.plot(ob[:, 0], ob[:, 1], color='#c8c8c8', lw=1.2, label='Sirius B (white dwarf, 1.02 suns)')
axo.plot(oa[:, 0], oa[:, 1], color='#7fb2ff', lw=2.2, label='Sirius A (2.06 suns)')
axo.scatter([0], [0], marker='+', s=120, color='#ffe24a', zorder=5, label='centre of mass')
for t, lab, col, ty in ((1844.2, '1844 and 1994: closest, 8 AU apart', '#ffd27a', None), (1862, '1862: Clark first sees B', '#ffd27a', -9.4), (1869.3, '1869 and 2019: farthest, 31.5 AU apart', '#bbb', -4.0), (2026, '2026', 'white', 1.2)):
    va, vb = top(posA(t) - bary(t)), top(posB(t) - bary(t))
    axo.plot([va[0], vb[0]], [va[1], vb[1]], color=col, lw=.8, ls=':')
    axo.scatter([va[0]], [va[1]], s=55, color=col, zorder=4); axo.scatter([vb[0]], [vb[1]], s=25, color=col, zorder=4)
    axo.annotate(lab, vb, xytext=(3.6, ty) if ty is not None else (-11.5, 3.4), textcoords='data', color=col, fontsize=10.5, va='center', arrowprops=dict(arrowstyle='-', color=col, lw=.6, alpha=.6, shrinkB=4))
axo.text(-11.5, -12.5, 'Each dotted line joins Sirius A (big dot) to Sirius B (small dot) at that date.\nSmall grey dots on B’s orbit: every 5 years.', color='#9aa', fontsize=9)
for t in range(1830, 1880, 5):      # little year ticks on B's orbit, 5 years apart
    vb = top(posB(t) - bary(t)); axo.scatter([vb[0]], [vb[1]], s=8, color='#888', zorder=3)
axo.set_xlim(-12, 14); axo.set_ylim(-13.5, 11); axo.set_aspect('equal')
axo.set_xlabel('AU toward galactic longitude 90° (Cygnus)  →', color='#ccc'); axo.set_ylabel('AU  (galactic centre is down)', color='#ccc')
axo.legend(loc='upper left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=10)
sun = top(-bary(1855)); sun = sun/np.linalg.norm(sun); mv = top(T @ vel); mv = mv/np.linalg.norm(mv)
cx, cy, L = 9.0, 8.0, 3.2
for vec, col, txt, off, ha, va in ((sun, '#ffe24a', 'toward the Sun (8.6 ly)', (0.3, 0), 'left', 'center'), (np.array([0, -1]), '#ff9f6b', 'galactic\ncentre', (0.25, -0.2), 'left', 'top'),
                                   (mv, '#7fb2ff', 'direction\nof travel', (-0.3, -0.2), 'right', 'top')):
    axo.annotate('', xy=(cx + L*vec[0], cy + L*vec[1]), xytext=(cx, cy), arrowprops=dict(arrowstyle='-|>', color=col, lw=1.8))
    axo.text(cx + L*vec[0] + off[0], cy + L*vec[1] + off[1], txt, color=col, fontsize=9.5, ha=ha, va=va)
axo.text(cx, cy + 0.5, 'directions', color='#9aa', fontsize=9, ha='center', va='bottom')

# 2. the whole 500 years, same view, AU from the Sun-ward end... relative to the 2026 centre of mass
C2 = top(bary(2026)); y5 = np.arange(1526, 2026.001, 0.2)
wa = np.array([top(posA(t)) - C2 for t in y5]); wb = np.array([top(bary(t)) - C2 for t in y5])
style(ax2, 'The whole 500 years, same view')
ax2.plot(wb[:, 0], wb[:, 1], color='#888', lw=1, ls=(0, (4, 3))); ax2.plot(wa[:, 0], wa[:, 1], color='#7fb2ff', lw=1.4)
for t in range(1550, 2026, 50):
    v = top(posA(t)) - C2; ax2.scatter([v[0]], [v[1]], s=26, color='#e6ecff', zorder=3)
    ax2.annotate(str(t), v, xytext=(8, 0), textcoords='offset points', color='#9fc4ff', fontsize=9, va='center')
v = top(posA(2026)) - C2; ax2.scatter([v[0]], [v[1]], s=160, color='white', zorder=4); ax2.annotate('2026', v, xytext=(9, -8), textcoords='offset points', color='white', fontsize=10.5, weight='bold')
x0, y0 = (top(posA(1830)) - C2); x1, y1 = (top(posA(1880)) - C2)
ax2.add_patch(plt.Rectangle((min(x0, x1) - 20, min(y0, y1) - 20), abs(x1 - x0) + 40, abs(y1 - y0) + 40, fill=False, ec='#c7a6ff', lw=1))
ax2.text(min(x0, x1) - 25, (y0 + y1)/2, 'shown at\nleft', color='#c7a6ff', fontsize=9, va='center', ha='right')
ax2.set_aspect('equal'); ax2.margins(0.15)
ax2.set_xlabel('AU toward l = 90°  →', color='#ccc'); ax2.set_ylabel('AU (galactic centre down)', color='#ccc')

fig.suptitle('Sirius A and B seen from the north galactic pole: their real path through space relative to the Sun', color='white', fontsize=15)
fig.text(.04, .015, 'Top view, as in the other galactic maps: looking down on the Milky Way’s disc from the north, Cygnus (l = 90°) to the right, the galactic centre down. '
         'Sirius moves at %.1f km/s relative to the Sun (%.1f km/s of it along the plane, %.1f km/s %s), so the pair advances about %d AU in each 50-year orbit while the stars circle '
         'about 20 AU apart. 1 AU = the Earth–Sun distance. The orbit is tilted 136° to our line of sight, so from the galactic pole it looks different from the sky view. '
         'Orbit: Bond et al. 2017; motion: Hipparcos (2007) and RV −5.5 km/s. The tilt’s front/back sense (which node is nearer us) is taken from the usual radial-velocity convention.'
         % (np.linalg.norm(vel)/KMS_LY, vt, abs(g[2])/KMS_LY, 'toward the south' if g[2] < 0 else 'toward the north', round(vt*KMS_LY*50*AU_LY, -1)),
         color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, 'sirius_motion_from_north_galactic_pole.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print('done')
