# The Milky Way seen from above (face-on), with galactic longitude drawn from the Sun and landmark stars and nebulae.
# Spiral arms: logarithmic spirals with the pitch angles and kinks fitted by Reid et al. (2019, ApJ 885, 131) to maser
# parallaxes, extended beyond the measured parts to give a whole-Galaxy picture: a schematic, not a photograph.
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
KPC = 3261.56                 # light years per kiloparsec
R0 = 8.15                     # Sun to galactic centre, kpc (26,600 ly)
# arm: (beta_kink deg, R_kink kpc, pitch below kink, pitch above kink, beta range drawn, colour, width kpc)
ARMS = {'Norma':           (18, 4.46, -1.0, 19.5, (-20, 330), '#8fb4ff', .30),
        'Scutum–Centaurus':      (23, 4.91, 14.1, 12.1, (-10, 340), '#a9c6ff', .42),
        'Sagittarius–Carina':    (24, 6.04, 17.1, 1.0, (-60, 200), '#9ec0ff', .30),
        'Perseus':               (40, 8.87, 10.3, 8.7, (-90, 240), '#b8d0ff', .42),
        'Local (Orion) Arm':     (9, 8.26, 11.4, 11.4, (-12, 34), '#ffe7a8', .22),
        'Outer':                 (18, 12.24, 3.0, 9.4, (-30, 90), '#7d9fe0', .32)}

def arm_xy(bk, rk, p1, p2, b0, b1, n=600):
    beta = np.linspace(b0, b1, n)
    psi = np.where(beta < bk, p1, p2)
    R = rk*np.exp(-np.radians(beta - bk)*np.tan(np.radians(psi)))
    # beta = 0 toward the Sun, increasing in the direction of galactic rotation (clockwise seen from the north)
    return R*np.sin(np.radians(beta)), R*np.cos(np.radians(beta)), R

def from_sun(l, d_ly):
    """galactocentric x, y (kpc) of an object at galactic longitude l and distance d from the Sun (in the plane)"""
    d = d_ly / KPC; lr = math.radians(l)
    return d*math.sin(lr), R0 - d*math.cos(lr)       # l = 0 points from the Sun to the centre; l = 90 toward rotation (+x)

# landmarks: name, galactic longitude, distance (ly)
FAR = [('Galactic centre (Sagittarius A*)', 0.0, 26600), ('Lagoon Nebula M8', 6.0, 4100), ('Eagle Nebula M16', 17.0, 5700),
       ('Cygnus X-1', 71.3, 7200), ('Deneb', 84.3, 2600), ('Cassiopeia A', 111.7, 11000), ('Double Cluster h & χ Persei', 134.6, 7500),
       ('Crab Nebula M1', 184.6, 6500), ('Orion Nebula M42', 209.0, 1340), ('Eta Carinae', 287.6, 7500), ('Omega Centauri', 309.1, 17100)]
NEAR = [('Hyades', 180.1, 153),
        ('Canopus', 261.2, 310), ('Polaris', 123.3, 433), ('Pleiades', 166.6, 444), ('Antares', 351.9, 550), ('Betelgeuse', 199.8, 550),
        ('Spica', 316.1, 250), ('Rigel', 209.2, 860), ('Orion Nebula', 209.0, 1340), ('Deneb', 84.3, 2600)]

fig = plt.figure(figsize=(16, 11), dpi=130); fig.patch.set_facecolor('#05070d')
ax = fig.add_axes([0.02, 0.06, 0.62, 0.86]); ax.set_facecolor('#05070d')
# soft disc glow: exponential disc + bulge, drawn as an image
g = np.linspace(-17, 17, 700); X, Y = np.meshgrid(g, g); RR = np.hypot(X, Y)
glow = np.exp(-RR/3.2) + 1.6*np.exp(-(RR/1.1)**2)
ax.imshow(glow, extent=(-17, 17, -17, 17), origin='lower', cmap='bone', alpha=.55, vmin=0, vmax=2.0, zorder=0)
# bar, about 27° from the Sun–centre line, half-length 4.5 kpc
ba = math.radians(27); bx = 4.5*math.sin(ba); by = 4.5*math.cos(ba)
ax.plot([-bx, bx], [-by, by], color='#ffe3b0', lw=14, alpha=.35, solid_capstyle='round', zorder=1)
for name, (bk, rk, p1, p2, (b0, b1), col, w) in ARMS.items():
    x, y, R = arm_xy(bk, rk, p1, p2, b0, b1)
    keep = R < 16.5
    ax.plot(x[keep], y[keep], color=col, lw=w*22, alpha=.28, solid_capstyle='round', zorder=2)
    ax.plot(x[keep], y[keep], color=col, lw=1.2, alpha=.9, zorder=3)
    bl = {'Perseus': -32, 'Outer': -22, 'Sagittarius–Carina': -38, 'Scutum–Centaurus': 200, 'Norma': 160, 'Local (Orion) Arm': 34}[name]
    lx_, ly_, _ = arm_xy(bk, rk, p1, p2, bl, bl, 1)
    ax.text(lx_[0], ly_[0], name, color=col, fontsize=9.5, style='italic', ha='center', va='bottom', zorder=6,
            bbox=dict(boxstyle='round,pad=.15', fc='#05070d', ec='none', alpha=.6))
# galactic longitude from the Sun
sx, sy = 0, R0
for l in range(0, 360, 30):
    lr = math.radians(l); L = 30 if l in (0,) else 40
    ex, ey = sx + L*math.sin(lr), sy - L*math.cos(lr)
    ax.plot([sx, ex], [sy, ey], color='#d9c27a', lw=.8 if l % 90 else 1.3, alpha=.75, zorder=4)
    # label where the line leaves the frame
    t = min([s for s in [(16.3 - sx)/math.sin(lr) if math.sin(lr) > 1e-9 else 1e9, (-16.3 - sx)/math.sin(lr) if math.sin(lr) < -1e-9 else 1e9,
                         (sy - (-16.3))/math.cos(lr) if math.cos(lr) > 1e-9 else 1e9, (sy - 16.3)/math.cos(lr) if math.cos(lr) < -1e-9 else 1e9] if s > 0])
    lx, ly = sx + t*math.sin(lr), sy - t*math.cos(lr)
    ax.text(lx, ly, 'l = %d°' % l, color='#f0d98a', fontsize=10 if l % 90 == 0 else 8.5, ha='center', va='center',
            bbox=dict(boxstyle='round,pad=.2', fc='#05070d', ec='none', alpha=.8), zorder=7)
for rr in (10000, 20000, 30000):
    ax.add_patch(Circle((sx, sy), rr/KPC, fill=False, ec='#d9c27a', lw=.6, ls=(0, (3, 4)), alpha=.6, zorder=4))
    ax.text(sx + rr/KPC*0.87, sy - rr/KPC*0.5, '%s ly' % format(rr, ','), color='#cdb874', fontsize=8, zorder=7)
for name, l, d in FAR:
    x, y = from_sun(l, d)
    ax.plot(x, y, marker='o', ms=5, color='#ff9d6e' if 'centre' not in name else '#ffffff', zorder=8)
    ax.annotate(name, (x, y), xytext=(6, 4), textcoords='offset points', color='#ffd2bd', fontsize=8.5, zorder=8)
ax.plot(sx, sy, marker='*', ms=16, color='#ffe24a', zorder=9); ax.annotate('Sun', (sx, sy), xytext=(8, -14), textcoords='offset points', color='#ffe24a', fontsize=11, weight='bold')
ax.annotate('', xy=(13.2, 7.0), xytext=(10.2, 11.6), arrowprops=dict(arrowstyle='->', color='#9fb6d8', lw=1.4, connectionstyle='arc3,rad=-.2'))
ax.text(12.9, 10.2, 'rotation\n(clockwise)', color='#9fb6d8', fontsize=9, ha='left')
ax.set_xlim(-16.5, 16.5); ax.set_ylim(-16.5, 16.5); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('The Milky Way seen from the north galactic pole, with galactic longitude measured from the Sun', color='white', fontsize=13)

# inset: the Sun's neighbourhood
ax2 = fig.add_axes([0.66, 0.30, 0.32, 0.46]); ax2.set_facecolor('#0b1020')
for s in ax2.spines.values(): s.set_color('#d9c27a')
M = 3000
for l in range(0, 360, 30):
    lr = math.radians(l); ax2.plot([0, M*1.3*math.sin(lr)], [0, -M*1.3*math.cos(lr)], color='#d9c27a', lw=.6, alpha=.6)
    ax2.text(M*.93*math.sin(lr), -M*.93*math.cos(lr), '%d°' % l, color='#f0d98a', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='round,pad=.15', fc='#0b1020', ec='none'))
for rr in (500, 1000, 2000):
    ax2.add_patch(Circle((0, 0), rr, fill=False, ec='#d9c27a', lw=.5, ls=(0, (3, 4)), alpha=.6)); ax2.text(rr*.71, rr*.71, '%s ly' % format(rr, ','), color='#cdb874', fontsize=7)
for name, l, d in NEAR:
    lr = math.radians(l); x, y = d*math.sin(lr), -d*math.cos(lr)
    ax2.plot(x, y, 'o', ms=4, color='#ff9d6e'); ax2.annotate(name, (x, y), xytext=(4, 3), textcoords='offset points', color='#ffd2bd', fontsize=7.5)
ax2.plot(0, 0, '*', ms=12, color='#ffe24a'); ax2.annotate('Sun (with Sirius, Vega, Arcturus\nand the Big Dipper core within 80 ly)', (0, 0), xytext=(-150, -30), textcoords='offset points', color='#ffe24a', fontsize=7.5)
ax2.set_xlim(-M, M); ax2.set_ylim(-M, M); ax2.set_aspect('equal'); ax2.set_xticks([]); ax2.set_yticks([])
ax2.set_title('Zoom: 3,000 light years around the Sun (same orientation)', color='white', fontsize=10)
ax2.text(0, -M*1.1, 'Galactic centre is straight down (l = 0°) in both views', color='#aaa', fontsize=8, ha='center')

fig.text(.02, .025, 'Galactic longitude l is measured in the plane from the Sun: 0° toward the galactic centre, 90° in the direction of the Sun’s orbit (Cygnus), '
         '180° away from the centre (Auriga), 270° toward Vela. Spiral arms: schematic after Reid et al. 2019 (maser parallaxes), extended beyond the measured parts; '
         'Sun 26,600 ly from the centre. Distances to landmarks are approximate.', color='#9aa', fontsize=8.5, wrap=True)
fig.savefig(os.path.join(here, 'milky_way_galactic_longitude_map.png'), facecolor=fig.get_facecolor(), bbox_inches='tight')
print('done')
