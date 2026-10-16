# The Pleiades and the Ursa Major nucleus close up, as seen from the direction of the galactic core.
# The viewer looks outward from the centre of the Galaxy (from l = 0° toward l = 180°): north galactic pole up,
# l = 270° (Vela) to the right and l = 90° (Cygnus) to the left; the mirror image of the view from the anticentre.
import csv, math, os, statistics
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
              [ 0.4941094279, -0.4448296300,  0.7469822445],
              [-0.8676661490, -0.1980763734,  0.4559837762]])

def xyz_radec(ra, dec, plx):
    a, d = math.radians(ra), math.radians(dec)
    return T @ (3261.56/plx*np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]))
def xyz_lb(l, b, ly):
    lr, br = math.radians(l), math.radians(b); return np.array([ly*math.cos(br)*math.cos(lr), ly*math.cos(br)*math.sin(lr), ly*math.sin(br)])

def colour_bprp(c):
    stops = [(-0.2, (0.62, 0.71, 1.0)), (0.3, (0.92, 0.94, 1.0)), (0.8, (1.0, 0.93, 0.62)), (1.4, (1.0, 0.70, 0.36)), (2.6, (1.0, 0.42, 0.30))]
    c = max(stops[0][0], min(stops[-1][0], c))
    for (x0, c0), (x1, c1) in zip(stops, stops[1:]):
        if c <= x1:
            t = (c - x0)/(x1 - x0); return tuple(c0[i] + t*(c1[i] - c0[i]) for i in range(3))
SPEC = {'O': -0.3, 'B': -0.1, 'A': 0.1, 'F': 0.5, 'G': 0.8, 'K': 1.2, 'M': 1.8}

# ---- Pleiades (Gaia members + the seven brightest, as in make_pleiades_views.py)
ple = []
for r in csv.DictReader(open(os.path.join(here, 'pleiades_members_cantat_gaudin_2020.csv'), encoding='utf-8')):
    try: ple.append(dict(ra=float(r['RA_ICRS']), dec=float(r['DE_ICRS']), plx=float(r['Plx']), m=float(r['Gmag']),
                         c=float(r['BP-RP']) if r['BP-RP'].strip() else 0.8, name=''))
    except ValueError: pass
plx0 = statistics.median(s['plx'] for s in ple)
for n, (ra, de) in {'Celaeno': (56.201, 24.289), 'Asterope': (56.477, 24.555)}.items():
    min(ple, key=lambda s: (s['ra'] - ra)**2 + (s['dec'] - de)**2)['name'] = n
for n, ra, de, V in [('Alcyone', 56.871, 24.105, 2.87), ('Atlas', 57.290, 24.053, 3.62), ('Electra', 56.219, 24.113, 3.70), ('Maia', 56.457, 24.368, 3.87),
                     ('Merope', 56.582, 23.948, 4.17), ('Taygeta', 56.302, 24.467, 4.29), ('Pleione', 57.297, 24.137, 5.05)]:
    ple.append(dict(ra=ra, dec=de, plx=plx0, m=V, c=-0.09, name=n))
for s in ple: s['p'] = xyz_radec(s['ra'], s['dec'], s['plx'])

# ---- Ursa Major nucleus and its neighbours (naked-eye stars, Hipparcos)
uma = []
for r in csv.DictReader(open(os.path.join(here, 'naked_eye_stars_within_250ly.csv'), encoding='utf-8')):
    uma.append(dict(p=xyz_lb(float(r['Galactic l']), float(r['Galactic b']), float(r['Distance (ly)'])), m=float(r['V mag']),
                    c=SPEC.get(('G' if r['Proper name'] == 'Capella' else r['Spectral type'])[:1].upper(), 0.6),
                    name=r['Proper name'] or (r['Name'] if r['Ursa Major group nucleus'] == 'yes' else ''), core=r['Ursa Major group nucleus'] == 'yes'))
core = [s for s in uma if s['core']]
cu = np.mean([s['p'] for s in core], axis=0)

def size(m, k): return max(2, min(300, (k - m)**2.2*0.9)) if m < k else 2

def panel(ax, rows, centre, R, kmag, title, rings, label_names, extra=None):
    ax.set_facecolor('#05070d')
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.tick_params(colors='#9aa', labelsize=8.5); ax.grid(color='#1c2636', lw=.6)
    # seen from the core: right = -Y (toward l = 270°), up = +Z; depth = X (nearer to the viewer = larger X): draw far ones first
    pts = [(s, s['p'] - centre) for s in rows]
    pts = [(s, d) for s, d in pts if abs(d[1]) < R and abs(d[2]) < R and abs(d[0]) < R]
    for s, d in sorted(pts, key=lambda t: t[1][0]):
        ax.scatter([-d[1]], [d[2]], s=size(s['m'], kmag), color=colour_bprp(s['c']), edgecolors='none', alpha=.95, zorder=3)
        if s['name'] and s['name'] in label_names:
            off = {'78 UMa': (6, -13), 'Alcor': (-6, -13), 'Alioth': (6, 5)}.get(s['name'], (5, 4))
            ax.annotate(s['name'], (-d[1], d[2]), xytext=off, textcoords='offset points', color='#f2f2f2', fontsize=9, zorder=5, ha='right' if off[0] < 0 else 'left')
    if extra: extra(ax, pts)
    for rr in rings: ax.add_patch(Circle((0, 0), rr, fill=False, ec='#9fb6d8', lw=.7, ls=(0, (3, 4)), alpha=.6))
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect('equal')
    ax.set_xlabel('light years toward l = 270° (Vela)  →', color='#ccc'); ax.set_ylabel('light years north of the galactic plane', color='#ccc')
    ax.set_title(title, color='white', fontsize=12.5)
    return len(pts)

fig = plt.figure(figsize=(19, 9.6), dpi=130); fig.patch.set_facecolor('#05070d')
a1 = fig.add_axes([0.04, 0.1, 0.43, 0.8]); a2 = fig.add_axes([0.54, 0.1, 0.43, 0.8])
cp = np.array([statistics.median(s['p'][k] for s in ple) for k in range(3)])
n1 = panel(a1, ple, cp, 45, 13.5, 'The Pleiades (%d stars), about 444 ly from the Sun' % len(ple), [38],
           {'Alcyone', 'Atlas', 'Electra', 'Maia', 'Merope', 'Taygeta', 'Celaeno', 'Pleione', 'Asterope'})
a1.text(26, -40, 'tidal radius ≈ 38 ly', color='#9fb6d8', fontsize=8.5)
def ring_core(ax, pts):
    for s, d in pts:
        if s['core']: ax.add_patch(Circle((-d[1], d[2]), 1.3, fill=False, ec='#d0302a', lw=1.2, zorder=4))
n2 = panel(a2, uma, cu, 28, 9.0, 'The Ursa Major nucleus and its neighbours, about 80 ly from the Sun',
           [], {s['name'] for s in uma if s['name']}, ring_core)
a2.text(-27, -26.5, 'red rings: the seven nucleus stars', color='#d0302a', fontsize=9)
for ax, c in ((a1, cp), (a2, cu)):
    # direction to the Sun projected on this view
    v = -c; dx, dy = -v[1], v[2]; n = math.hypot(dx, dy); R = ax.get_xlim()[1]
    ax.annotate('', xy=(dx/n*R*.93, dy/n*R*.93), xytext=(dx/n*R*.68, dy/n*R*.68), arrowprops=dict(arrowstyle='-|>', color='#ffe24a', lw=1.8))
    ax.text(dx/n*R*.62, dy/n*R*.62, 'Sun’s direction', color='#ffe24a', fontsize=9, ha='center')
fig.suptitle('Seen from the galactic core, zoomed in: looking outward from the centre of the Milky Way (north galactic pole up)', color='white', fontsize=14.5)
fig.text(.04, .02, 'An orthographic close-up in the direction from the galactic centre (26,600 ly away) toward each group. From the core both groups would appear as faint specks; '
         'this is the same direction of view, zoomed in. It is the mirror image of the edge-on view from the anticentre: here Vela (l = 270°) is to the right and Cygnus (l = 90°) to the left. '
         'Pleiades: Gaia members (Cantat-Gaudin & Anders 2020); Ursa Major: naked-eye Hipparcos stars within 28 ly of the nucleus centre.', color='#9aa', fontsize=9, wrap=True)
fig.savefig(os.path.join(here, 'pleiades_and_ursa_major_from_galactic_core.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)
print('done', n1, n2)
