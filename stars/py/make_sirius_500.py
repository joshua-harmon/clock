# Sirius's proper motion over the last 500 years (1526 to 2026), zoomed in: offsets on the sky in arcminutes,
# fixed J2000 coordinates (no precession), straight-line motion from Hipparcos (2007).
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
K = 4.740470446; KMS_LY = 1.0227e-6*3.26156
RA0, DE0, PLX, PMRA, PMDE, RV = 101.2872, -16.7161, 379.21, -546.01, -1223.07, -5.5   # J2000 position, epoch 1991.25
EP0 = 1991.25

def unit(ra, de):
    a, d = math.radians(ra), math.radians(de); return np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
r0 = unit(RA0, DE0); a, d = math.radians(RA0), math.radians(DE0)
p = np.array([-math.sin(a), math.cos(a), 0]); q = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
D0 = 3261.56/PLX
vel = ((K*PMRA/PLX)*p + (K*PMDE/PLX)*q + RV*r0)*KMS_LY

def offset(year):
    """arcminutes east and north of the star's 2026 place"""
    def pos(y):
        v = D0*r0 + (y - EP0)*vel; return v/np.linalg.norm(v)
    u, u26 = pos(year), pos(2026)
    return 60*math.degrees(np.dot(u - u26, -p)) * -1, 60*math.degrees(np.dot(u - u26, q))   # east positive

yrs = np.arange(1526, 2027, 1)
xy = np.array([offset(y) for y in yrs])
tot = math.hypot(*xy[0])
print('moved %.2f arcmin in 500 years' % tot)

fig, ax = plt.subplots(figsize=(11, 11), dpi=130); fig.patch.set_facecolor('#05070d'); ax.set_facecolor('#05070d')
for s in ax.spines.values(): s.set_color('#3a4a60')
ax.tick_params(colors='#9aa'); ax.grid(color='#1c2636', lw=.6)
ax.plot(xy[:, 0], xy[:, 1], color='#7fb2ff', lw=1.6, zorder=2)
for y in range(1526, 2027, 50):
    x_, y_ = offset(y); ax.scatter([x_], [y_], s=40, color='#e6ecff', edgecolors='#7fb2ff', zorder=3)
    if y != 2026: ax.annotate(str(y), (x_, y_), xytext=(9, -3), textcoords='offset points', color='#9fc4ff', fontsize=9.5)
for y, txt, dy in [(1718, '1718: Edmond Halley discovers that Sirius, Arcturus\nand Aldebaran have moved since Ptolemy’s catalogue:\nthe first proper motions', (-6, 13)),
               (1844, '1844: Friedrich Bessel finds Sirius wobbling along its path\nand infers an unseen companion', (-6, 7.5)),
               (1862, '1862: Alvan Clark first sees the companion, Sirius B', (-6, 3.5))]:
    x_, y_ = offset(y); ax.scatter([x_], [y_], s=55, color='#ffd27a', zorder=4)
    ax.annotate(txt, (x_, y_), xytext=dy, textcoords='data', color='#ffd27a', fontsize=9.5, ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color='#ffd27a', lw=.7, alpha=.7, shrinkB=4))
x26, y26 = offset(2026)
ax.scatter([x26], [y26], s=220, color='white', zorder=5); ax.annotate('Sirius in 2026', (x26, y26), xytext=(12, -5), textcoords='offset points', color='white', fontsize=12, weight='bold')
# scale: the full Moon (31 arcmin across)
ax.add_patch(Circle((-14, -16), 15.5, fill=False, ec='#888', lw=1, ls=(0, (4, 3))))
ax.text(-14, -16, 'size of the\nfull Moon', color='#999', fontsize=9, ha='center', va='center')
ax.set_xlim(20, -32); ax.set_ylim(-34, 18); ax.set_aspect('equal')    # east to the left, as on a sky map
ax.set_xlabel('arcminutes east of today’s position (east to the left)', color='#ccc'); ax.set_ylabel('arcminutes north of today’s position', color='#ccc')
ax.set_title('Sirius’s proper motion, 1526 to 2026: %.1f arcminutes, about a third of the full Moon’s width' % tot, color='white', fontsize=13)
fig.text(.06, .015, 'Hipparcos (2007) proper motion, 1.339″ a year toward position angle 204° (south-southwest), with a radial velocity of −5.5 km/s; fixed J2000 coordinates, '
         'so precession of the equinoxes is not included. Sirius A’s 50-year orbit around Sirius B adds a wobble of a few arcseconds, too small to see at this scale.',
         color='#9aa', fontsize=8.5, wrap=True)
fig.savefig(os.path.join(here, 'sirius_motion_500_years.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)
print('done')
