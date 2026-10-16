# The Pleiades zoomed in on the named stars, seen three ways (from Earth, from the north galactic pole, edge-on).
# Same layout as make_pleiades_views.py, but zoomed so that only the named stars matter; text stays the same size.
# Members: Lodieu et al. 2019 (A&A 628, A66), Gaia DR2, members within the tidal radius (stars/pleiades/pleiades_lodieu2019_gaiadr2.tsv).
# The named stars are too bright for reliable Gaia DR2 parallaxes, so they are placed at the cluster's median distance.
import math, os, statistics
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
tsv = os.path.join(here, '..', 'pleiades', 'pleiades_lodieu2019_gaiadr2.tsv'); out = os.path.join(here, '..', '..', 'images', 'U6')
T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]])
NAMED = [('Alcyone', 56.871, 24.105, 2.87), ('Atlas', 57.290, 24.053, 3.62), ('Electra', 56.219, 24.113, 3.70), ('Maia', 56.457, 24.368, 3.87),
         ('Merope', 56.582, 23.948, 4.17), ('Taygeta', 56.302, 24.467, 4.29), ('Pleione', 57.297, 24.137, 5.05), ('Celaeno', 56.201, 24.289, 5.45),
         ('Asterope', 56.477, 24.555, 5.76)]
rows = [l.rstrip('\n').split('\t') for l in open(tsv) if l.strip() and not l.startswith('#')]
H = rows[0]; ix = {k: H.index(k) for k in H}
stars = []
for r in rows[1:]:
    if len(r) != len(H) or not r[ix['RA_ICRS']].strip().replace('.', '').isdigit(): continue
    if float(r[ix['distC']]) > 11.6: continue
    stars.append(dict(ra=float(r[ix['RA_ICRS']]), dec=float(r[ix['DE_ICRS']]), plx=float(r[ix['PlxCP']]), G=float(r[ix['GmagCP']] or 20), name=''))
plx0 = statistics.median(s['plx'] for s in stars)
named = [dict(ra=a, dec=d, plx=plx0, G=V, name=n) for n, a, d, V in NAMED]
def xyz(s):
    a, d = math.radians(s['ra']), math.radians(s['dec']); dist = 3261.56/s['plx']
    return T @ (dist*np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]))
for s in stars + named: s['X'], s['Y'], s['Z'] = xyz(s)
C = np.array([np.mean([s[k] for s in named]) for k in 'XYZ']); D0 = np.linalg.norm(C)
ra0 = np.mean([s['ra'] for s in named]); de0 = np.mean([s['dec'] for s in named])
def size(G): return max(1.5, min(420, (11 - G)**2.4*1.6)) if G < 11 else 1.5
BG = '#05070d'
fig = plt.figure(figsize=(21, 8.6), dpi=130); fig.patch.set_facecolor(BG)
panels = [fig.add_axes([0.01 + k*0.33, 0.1, 0.31, 0.78]) for k in range(3)]
for ax in panels:
    ax.set_facecolor(BG)
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.tick_params(colors='#9aa', labelsize=9); ax.grid(color='#1c2636', lw=.6)
OFF = {}
def draw(ax, fx, fy, lim):
    for s in sorted(stars, key=lambda s: -s['G']):
        x, y = fx(s), fy(s)
        if abs(x) < lim and abs(y) < lim: ax.scatter([x], [y], s=3, color='#5d6a85', edgecolors='none')
    for s in named:
        x, y = fx(s), fy(s); ax.scatter([x], [y], s=size(s['G']), color='#e8f0ff', edgecolors='#9fc0ff', linewidths=.6, zorder=5)
        o = OFF.get(s['name'], (9, 7))
        ax.annotate(s['name'], (x, y), xytext=o, textcoords='offset points', color='#ffd98a', fontsize=11, weight='bold', zorder=6, ha='right' if o[0] < 0 else 'left')
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
# 1. sky
L1 = 0.75
OFF = {}
draw(panels[0], lambda s: -(s['ra'] - ra0)*math.cos(math.radians(de0)), lambda s: s['dec'] - de0, L1)
panels[0].set_xlabel('degrees east ← → west', color='#ccc', fontsize=10); panels[0].set_ylabel('degrees north', color='#ccc', fontsize=10)
panels[0].set_title('From Earth: the named stars on the sky (1.5° across)', color='white', fontsize=12)
# 2, 3. 3-D views centred on the named stars
L2 = 6
for ax, (hk, vk, sgn, xl, yl, ttl) in zip(panels[1:], [
        ('Y', 'X', -1, 'light years toward l = 90° (Cygnus)', 'light years away from the galactic centre', 'From above: seen from the north galactic pole (12 ly across)'),
        ('Y', 'Z', 1, 'light years toward l = 90° (Cygnus)', 'light years north of the galactic plane', 'Turned 90°: edge-on from the galactic anticentre (12 ly across)')]):
    ih, iv = 'XYZ'.index(hk), 'XYZ'.index(vk)
    OFF = {'Atlas': (-10, -14), 'Pleione': (9, 7), 'Maia': (-10, 7), 'Asterope': (9, 8), 'Taygeta': (9, -2), 'Celaeno': (9, -12), 'Electra': (-10, -14)}
    draw(ax, lambda s, ih=ih: s[hk] - C[ih], lambda s, iv=iv: sgn*(s[vk] - C[iv]), L2)
    dx, dy = -C[ih], sgn*(-C[iv]); n = math.hypot(dx, dy)
    ax.annotate('', xy=(dx/n*L2*.92, dy/n*L2*.92), xytext=(dx/n*L2*.62, dy/n*L2*.62), arrowprops=dict(arrowstyle='-|>', color='#ffe24a', lw=1.8))
    ax.text(dx/n*L2*.5 + L2*.12, dy/n*L2*.5, 'toward the Sun\n(%d ly)' % round(D0, -1), color='#ffe24a', fontsize=9.5, ha='center')
    ax.set_xlabel(xl, color='#ccc', fontsize=10); ax.set_ylabel(yl, color='#ccc', fontsize=10); ax.set_title(ttl, color='white', fontsize=12)
fig.suptitle('The named stars of the Pleiades, about %d light years away, seen three ways' % round(D0, -1), color='white', fontsize=15)
fig.text(.01, .02, 'Named stars in white, labelled; small grey dots are other Gaia members in the same field (Lodieu et al. 2019, Gaia DR2). The named stars are too bright for '
         'reliable Gaia parallaxes, so all nine are placed at the cluster’s median distance; the two 3-D views therefore show their spread across the sky but not '
         'their true depth along the line toward the Sun.', color='#9aa', fontsize=9, wrap=True)
fig.savefig(os.path.join(out, 'pleiades_named_stars_three_views.png'), facecolor=BG, bbox_inches='tight'); print('saved', len(stars))
