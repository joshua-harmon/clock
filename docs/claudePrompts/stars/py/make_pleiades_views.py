# The Pleiades (M45) star by star: as seen from Earth, from the north galactic pole, and edge-on.
# Members: Cantat-Gaudin & Anders (2020, A&A 640, A1), Gaia DR2, "Melotte 22" (952 stars, G 5.4 to 18).
# The seven brightest stars are too bright for that list; they are added at their sky positions and placed
# at the cluster's average distance (their own parallaxes are less reliable).
import csv, math, os, statistics
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
              [ 0.4941094279, -0.4448296300,  0.7469822445],
              [-0.8676661490, -0.1980763734,  0.4559837762]])
BRIGHT = [('Alcyone', 56.871, 24.105, 2.87, -0.09), ('Atlas', 57.290, 24.053, 3.62, -0.08), ('Electra', 56.219, 24.113, 3.70, -0.11),
          ('Maia', 56.457, 24.368, 3.87, -0.07), ('Merope', 56.582, 23.948, 4.17, -0.06), ('Taygeta', 56.302, 24.467, 4.29, -0.11),
          ('Pleione', 57.297, 24.137, 5.05, -0.08)]
NAMED = {'Celaeno': (56.201, 24.289), 'Asterope': (56.477, 24.555)}

stars = []
for r in csv.DictReader(open(os.path.join(here, 'pleiades_members_cantat_gaudin_2020.csv'), encoding='utf-8')):
    try: stars.append(dict(ra=float(r['RA_ICRS']), dec=float(r['DE_ICRS']), plx=float(r['Plx']), G=float(r['Gmag']),
                           c=float(r['BP-RP']) if r['BP-RP'].strip() else 0.8, name='', assumed=False))
    except ValueError: pass
plx0 = statistics.median(s['plx'] for s in stars)
for n, (ra, de) in NAMED.items():
    min(stars, key=lambda s: (s['ra'] - ra)**2 + (s['dec'] - de)**2)['name'] = n
for n, ra, de, V, c in BRIGHT:
    stars.append(dict(ra=ra, dec=de, plx=plx0, G=V, c=c, name=n, assumed=True))

for s in stars:
    a, d = math.radians(s['ra']), math.radians(s['dec']); dist = 3261.56/s['plx']
    s['X'], s['Y'], s['Z'] = T @ (dist*np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]))
    s['ly'] = dist
C = np.array([statistics.median(s[k] for s in stars) for k in 'XYZ'])
D0 = np.linalg.norm(C)
l0 = math.degrees(math.atan2(C[1], C[0])) % 360; b0 = math.degrees(math.asin(C[2]/D0))
print('%d stars; centre %.0f ly at l = %.1f, b = %.1f' % (len(stars), D0, l0, b0))

def colour(c):
    stops = [(-0.2, (0.62, 0.71, 1.0)), (0.3, (0.92, 0.94, 1.0)), (0.8, (1.0, 0.93, 0.62)), (1.4, (1.0, 0.70, 0.36)), (2.6, (1.0, 0.42, 0.30))]
    c = max(stops[0][0], min(stops[-1][0], c))
    for (x0, c0), (x1, c1) in zip(stops, stops[1:]):
        if c <= x1:
            t = (c - x0)/(x1 - x0); return tuple(c0[i] + t*(c1[i] - c0[i]) for i in range(3))
def size(G): return max(2, min(260, (13.5 - G)**2.2*0.9)) if G < 13.5 else 2

order = sorted(stars, key=lambda s: -s['G'])
fig = plt.figure(figsize=(21, 8.6), dpi=130); fig.patch.set_facecolor('#05070d')
panels = [fig.add_axes([0.01 + k*0.33, 0.1, 0.31, 0.78]) for k in range(3)]
for ax in panels:
    ax.set_facecolor('#05070d')
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.tick_params(colors='#9aa', labelsize=8.5); ax.grid(color='#1c2636', lw=.6)

# 1. from Earth: the sky, degrees from the cluster centre (east to the left, north up)
ax = panels[0]; ra0, de0 = statistics.median(s['ra'] for s in stars), statistics.median(s['dec'] for s in stars)
for s in order:
    x = -(s['ra'] - ra0)*math.cos(math.radians(de0)); y = s['dec'] - de0
    ax.scatter([x], [y], s=size(s['G']), color=colour(s['c']), edgecolors='none')
    if s['name']: ax.annotate(s['name'], (x, y), xytext=(5, 4), textcoords='offset points', color='#f2f2f2', fontsize=8.5)
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_aspect('equal')
ax.set_xlabel('degrees east ← → west', color='#ccc'); ax.set_ylabel('degrees north', color='#ccc')
ax.set_title('From Earth: the Pleiades on the sky (6° across)', color='white', fontsize=12)

# 2. from the north galactic pole, centred on the cluster; galactic centre down, as in the other maps
def sunarrow(ax, dx, dy, R):
    n = math.hypot(dx, dy); ax.annotate('', xy=(dx/n*R*.92, dy/n*R*.92), xytext=(dx/n*R*.62, dy/n*R*.62), arrowprops=dict(arrowstyle='-|>', color='#ffe24a', lw=1.8))
    ax.text(dx/n*R*.6, dy/n*R*.6 + R*.05, 'toward the Sun\n(%d ly)' % round(D0, -1), color='#ffe24a', fontsize=9, ha='center')
R = 45
for ax, (hk, vk, sgn, xl, yl, ttl) in zip(panels[1:], [
        ('Y', 'X', -1, 'light years toward l = 90° (Cygnus)', 'light years away from the galactic centre', 'From above: seen from the north galactic pole'),
        ('Y', 'Z', 1, 'light years toward l = 90° (Cygnus)', 'light years north of the galactic plane', 'Turned 90°: seen edge-on from the galactic anticentre')]):
    ih, iv = 'XYZ'.index(hk), 'XYZ'.index(vk)
    for s in order:
        x = s[hk] - C[ih]; y = sgn*(s[vk] - C[iv])
        ax.scatter([x], [y], s=size(s['G']), color=colour(s['c']), edgecolors='none', alpha=.95)
        if s['name'] in ('Alcyone', 'Atlas', 'Electra', 'Maia', 'Merope', 'Taygeta', 'Celaeno'):
            ax.annotate(s['name'], (x, y), xytext=(5, 4), textcoords='offset points', color='#f2f2f2', fontsize=8.5)
    ax.add_patch(Circle((0, 0), 38, fill=False, ec='#9fb6d8', lw=.8, ls=(0, (3, 4)), alpha=.7))
    ax.text(27, -36, 'tidal radius ≈ 38 ly', color='#9fb6d8', fontsize=8.5)
    sunarrow(ax, -C[ih], sgn*(-C[iv]), R)
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect('equal')
    ax.set_xlabel(xl, color='#ccc'); ax.set_ylabel(yl, color='#ccc'); ax.set_title(ttl, color='white', fontsize=12)
fig.suptitle('The Pleiades, %d stars, about %d light years away, seen three ways' % (len(stars), round(D0, -1)), color='white', fontsize=15)
fig.text(.01, .02, 'Members from Gaia (Cantat-Gaudin & Anders 2020); colour from the stars’ Gaia colour (blue-white hot stars to red dwarfs), dot size by brightness. '
         'The seven brightest stars are placed at the cluster’s average distance. In the two 3-D views, distance errors of a few light years stretch the cluster slightly '
         'along the line toward the Sun. The cluster is about 100 million years old and is slowly dissolving into the Galaxy.', color='#9aa', fontsize=9, wrap=True)
fig.savefig(os.path.join(here, 'pleiades_three_views.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)
print('done')
