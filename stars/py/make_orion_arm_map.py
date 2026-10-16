# The Orion (Local) Arm around the Sun, seen from the north galactic pole, with galactic longitude from the Sun
# and landmark stars, clusters and nebulae at their measured directions and approximate distances.
# Arms: logarithmic spirals after Reid et al. (2019, ApJ 885, 131); the band widths are the fitted arm widths.
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
KPC = 3261.56
R0 = 8.15
ARMS = {'Sagittarius–Carina Arm': (24, 6.04, 17.1, 1.0, 0.27, '#8fb4ff'),
        'Local (Orion) Arm':      (9, 8.26, 11.4, 11.4, 0.31, '#ffe7a8'),
        'Perseus Arm':            (40, 8.87, 10.3, 8.7, 0.35, '#b8d0ff')}
M = 8500     # half-width of the map, light years

def arm(bk, rk, p1, p2, b0=-80, b1=80, n=1500):
    beta = np.linspace(b0, b1, n); psi = np.where(beta < bk, p1, p2)
    R = rk*np.exp(-np.radians(beta - bk)*np.tan(np.radians(psi)))
    gx, gy = R*np.sin(np.radians(beta)), R*np.cos(np.radians(beta))
    return gx*KPC, (gy - R0)*KPC                     # Sun-centred: x toward l = 90°, y away from the centre (l = 180°)

def pos(l, d):
    lr = math.radians(l); return d*math.sin(lr), -d*math.cos(lr)

# name, galactic longitude (deg), distance (ly), kind
LAND = [('Betelgeuse', 199.8, 550, 'star'), ('Rigel', 209.2, 860, 'star'), ('Antares', 351.9, 550, 'star'), ('Deneb', 84.3, 2600, 'star'),
        ('Canopus', 261.2, 310, 'star'), ('Polaris', 123.3, 433, 'star'), ('Spica', 316.1, 250, 'star'), ('Alnilam', 205.2, 2000, 'star'),
        ('Sadr', 78.1, 1800, 'star'), ('P Cygni', 75.8, 5100, 'star'), ('Mu Cephei (Garnet Star)', 100.6, 3000, 'star'),
        ('VY Canis Majoris', 239.4, 3900, 'star'), ('Wezen', 238.4, 1600, 'star'), ('Naos (Zeta Puppis)', 256.0, 1080, 'star'),
        ('Gamma Velorum', 262.8, 1100, 'star'), ('Eta Carinae', 287.6, 7500, 'star'),
        ('Pleiades', 166.6, 444, 'cluster'), ('Hyades', 180.1, 153, 'cluster'), ('Beehive (Praesepe)', 205.9, 577, 'cluster'),
        ('Double Cluster', 134.6, 7500, 'cluster'), ('IC 2391', 270.4, 500, 'cluster'),
        ('Orion Nebula', 209.0, 1344, 'nebula'), ('Rosette Nebula', 206.3, 5200, 'nebula'), ('North America Nebula', 85.0, 2600, 'nebula'),
        ('Lagoon Nebula', 6.0, 4100, 'nebula'), ('Trifid Nebula', 7.0, 5200, 'nebula'), ('Omega Nebula M17', 15.1, 5500, 'nebula'),
        ('Eagle Nebula', 17.0, 5700, 'nebula'), ('Crab Nebula', 184.6, 6500, 'nebula'), ('Vela supernova remnant', 263.9, 800, 'nebula'),
        ('Coalsack', 301.0, 600, 'nebula'), ('Rho Ophiuchi cloud', 353.0, 450, 'nebula'), ('Taurus dark clouds', 172.0, 450, 'nebula'),
        ('Cygnus X-1', 71.3, 7200, 'nebula')]
STY = {'star': ('o', '#ffcf6e', 6), 'cluster': ('D', '#9fe0ff', 5.5), 'nebula': ('s', '#ff8fb0', 5)}
OFF = {'Rigel': (6, -10), 'Alnilam': (6, -2), 'Orion Nebula': (-8, 6), 'Betelgeuse': (6, 2), 'Beehive (Praesepe)': (-10, 8),
       'Trifid Nebula': (6, 2), 'Lagoon Nebula': (6, -10), 'Omega Nebula M17': (6, 2), 'Eagle Nebula': (6, -8), 'Wezen': (6, -10),
       'Gamma Velorum': (6, 4), 'Naos (Zeta Puppis)': (-8, -12), 'Sadr': (-8, 6), 'North America Nebula': (6, -12), 'Deneb': (6, 4),
       'Pleiades': (-8, 6), 'Taurus dark clouds': (-8, -12), 'Hyades': (6, -8), 'IC 2391': (6, 4), 'Vela supernova remnant': (6, -10),
       'Canopus': (-8, -12), 'Polaris': (6, 2), 'Spica': (6, -8), 'Antares': (6, 2), 'Rho Ophiuchi cloud': (-8, -12)}
LEFT = {'Orion Nebula', 'Beehive (Praesepe)', 'Naos (Zeta Puppis)', 'Sadr', 'Pleiades', 'Taurus dark clouds', 'Canopus', 'Rho Ophiuchi cloud'}

fig = plt.figure(figsize=(21, 14), dpi=130); fig.patch.set_facecolor('#05070d')
ax = fig.add_axes([0.01, 0.07, 0.63, 0.88]); ax.set_facecolor('#05070d')
for name, (bk, rk, p1, p2, w, col) in ARMS.items():
    x, y = arm(bk, rk, p1, p2)
    k = (abs(x) < 1.6*M) & (abs(y) < 1.6*M)
    ax.plot(x[k], y[k], color=col, lw=2*w*KPC/(2*M)*887, alpha=.10, solid_capstyle='butt', zorder=1)   # band about the fitted width
    ax.plot(x[k], y[k], color=col, lw=1.0, alpha=.7, ls=(0, (6, 4)), zorder=2)
    # label where the arm crosses the left part of the map
    i = np.argmin(abs(x - (-0.62*M)) + 1e6*(abs(y) > 0.95*M))
    ax.text(x[i], y[i] + 250, name, color=col, fontsize=12, style='italic', ha='center', va='bottom', zorder=6,
            bbox=dict(boxstyle='round,pad=.2', fc='#05070d', ec='none', alpha=.7))
for l in range(0, 360, 15):
    lr = math.radians(l); L = M*1.6
    major = l % 30 == 0
    ax.plot([0, L*math.sin(lr)], [0, -L*math.cos(lr)], color='#d9c27a', lw=1.1 if l % 90 == 0 else (.7 if major else .35), alpha=.75 if major else .4, zorder=3)
    if major:
        s, c = math.sin(lr), -math.cos(lr); t = min(M*0.965/abs(s) if abs(s) > 1e-9 else 1e9, M*0.965/abs(c) if abs(c) > 1e-9 else 1e9)
        ax.text(t*s, t*c, '%d°' % l, color='#f0d98a', fontsize=11 if l % 90 == 0 else 9.5, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=.2', fc='#05070d', ec='none', alpha=.85), zorder=7)
for rr in (1000, 2000, 4000, 6000, 8000):
    ax.add_patch(Circle((0, 0), rr, fill=False, ec='#d9c27a', lw=.6, ls=(0, (3, 4)), alpha=.55, zorder=3))
    ax.text(-rr*0.97, rr*0.26 + 60, '%s ly' % format(rr, ','), color='#cdb874', fontsize=8.5, zorder=7)
for name, l, d, kind in LAND:
    x, y = pos(l, d); mk, col, ms = STY[kind]
    ax.plot(x, y, marker=mk, ms=ms, color=col, zorder=8, ls='none')
    if d >= 1000:
        ox, oy = OFF.get(name, (6, 3))
        ax.annotate(name, (x, y), xytext=(ox, oy), textcoords='offset points', color=col, fontsize=9.5, zorder=9,
                    ha='right' if name in LEFT else 'left')
ax.plot(0, 0, marker='*', ms=17, color='#ffe24a', zorder=10)
ax.annotate('Sun', (0, 0), xytext=(9, -16), textcoords='offset points', color='#ffe24a', fontsize=12, weight='bold', zorder=10)
ax.text(0, -M*0.9, '↓  toward the galactic centre, 26,600 ly', color='#e8e8e8', fontsize=11, ha='center', zorder=7)
ax.text(0, M*0.88, 'away from the galactic centre  ↑', color='#e8e8e8', fontsize=11, ha='center', zorder=7)
ax.text(M*0.9, M*0.05, 'Sun’s orbit →\n(Cygnus)', color='#e8e8e8', fontsize=10, ha='right', zorder=7)
ax.text(-M*0.9, M*0.05, '(Vela)', color='#e8e8e8', fontsize=10, ha='left', zorder=7)
for k, (lab, kind) in enumerate([('star', 'star'), ('star cluster', 'cluster'), ('nebula, cloud or remnant', 'nebula')]):
    mk, col, ms = STY[kind]; ax.plot(-M*0.95, -M*(0.74 + k*0.045), marker=mk, ms=ms, color=col); ax.text(-M*0.92, -M*(0.74 + k*0.045), lab, color='#ccc', fontsize=9.5, va='center')
ax.set_xlim(-M, M); ax.set_ylim(-M, M); ax.set_aspect('equal'); ax.axis('off')
ax.add_patch(plt.Rectangle((-1250, -1250), 2500, 2500, fill=False, ec='#9fb6d8', lw=1, zorder=9))
ax.text(1300, -1250, 'zoom at right →', color='#9fb6d8', fontsize=9, va='bottom', zorder=9)
# inset: 1,250 ly around the Sun, same orientation
Z = 1250
ax2 = fig.add_axes([0.655, 0.24, 0.33, 0.495]); ax2.set_facecolor('#0b1020')
for sp in ax2.spines.values(): sp.set_color('#9fb6d8')
for nm, (bk, rk, p1, p2, w, col) in ARMS.items():
    x, y = arm(bk, rk, p1, p2, -15, 15, 600); ax2.plot(x, y, color=col, lw=1, ls=(0, (6, 4)), alpha=.6)
for l in range(0, 360, 30):
    lr = math.radians(l); ax2.plot([0, 2*Z*math.sin(lr)], [0, -2*Z*math.cos(lr)], color='#d9c27a', lw=.5, alpha=.5)
    s_, c_ = math.sin(lr), -math.cos(lr); t = min(Z*.93/abs(s_) if abs(s_) > 1e-9 else 1e9, Z*.93/abs(c_) if abs(c_) > 1e-9 else 1e9)
    ax2.text(t*s_, t*c_, '%d°' % l, color='#f0d98a', fontsize=7.5, ha='center', va='center', bbox=dict(boxstyle='round,pad=.12', fc='#0b1020', ec='none'))
for rr in (250, 500, 1000):
    ax2.add_patch(Circle((0, 0), rr, fill=False, ec='#d9c27a', lw=.5, ls=(0, (3, 4)), alpha=.5)); ax2.text(rr*.97, -rr*.26, '%d ly' % rr, color='#cdb874', fontsize=7)
ZO = {'Rigel': (5, 3), 'Orion Nebula': (-6, 4), 'Betelgeuse': (-6, -10), 'Beehive (Praesepe)': (-6, 3), 'Pleiades': (5, 3), 'Taurus dark clouds': (5, -10),
      'Hyades': (5, 2), 'Polaris': (5, 2), 'Canopus': (5, 4), 'IC 2391': (5, -10), 'Vela supernova remnant': (-6, 5), 'Coalsack': (-6, -9),
      'Spica': (5, -9), 'Antares': (5, -3), 'Rho Ophiuchi cloud': (-6, -9)}
for name, l, d, kind in LAND:
    if d > 1300: continue
    x, y = pos(l, d); mk, col, ms = STY[kind]
    ax2.plot(x, y, marker=mk, ms=ms - 1, color=col, ls='none')
    ox, oy = ZO.get(name, (5, 3))
    ax2.annotate(name, (x, y), xytext=(ox, oy), textcoords='offset points', color=col, fontsize=7.5, ha='right' if ox < 0 else 'left')
ax2.plot(0, 0, '*', ms=12, color='#ffe24a'); ax2.annotate('Sun', (0, 0), xytext=(5, -10), textcoords='offset points', color='#ffe24a', fontsize=8)
ax2.set_xlim(-Z, Z); ax2.set_ylim(-Z, Z); ax2.set_aspect('equal'); ax2.set_xticks([]); ax2.set_yticks([])
ax2.set_title('Zoom: 1,250 ly around the Sun', color='white', fontsize=10)
ax.set_title('The Orion Arm: 8,500 light years around the Sun, galactic longitude from the Sun', color='white', fontsize=15)
fig.text(.03, .03, 'Seen from the north galactic pole; positions are projected onto the galactic plane. Galactic longitude: 0° toward the galactic centre, 90° in the direction of the Sun’s orbit, '
         '180° away from the centre, 270° toward Vela. Dashed lines and shaded bands: spiral arms after Reid et al. 2019 (maser parallaxes). '
         'Distances are approximate, especially beyond 2,000 ly.', color='#9aa', fontsize=9, wrap=True)
fig.savefig(os.path.join(here, 'orion_arm_galactic_longitude_map.png'), facecolor=fig.get_facecolor(), bbox_inches='tight')
print('done')
