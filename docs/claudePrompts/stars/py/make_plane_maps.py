# Galactic-plane maps in the style of the Atlas of the Universe: the plane of the Milky Way seen at a slant,
# each star on a stalk rising from the plane (solid above, dashed below), rings of distance, and the directions
# in the plane marked by galactic longitude (0° = toward the centre of the Galaxy).
# Reads stars_within_10ly.csv, stars_within_50ly.csv and naked_eye_stars_within_250ly.csv (same folder).
import csv, math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

here = os.path.dirname(os.path.abspath(__file__))
ELEV = 32            # viewing angle above the plane, degrees
SE, CE = math.sin(math.radians(ELEV)), math.cos(math.radians(ELEV))
COL = {'O': '#9bb0ff', 'B': '#aabfff', 'A': '#e6ecff', 'F': '#fff8e8', 'G': '#ffe680', 'K': '#ffb056', 'M': '#ff6a4d', 'D': '#ffffff'}
DIRS = [(0, 'galactic centre\n(Sagittarius)'), (90, 'direction of the Sun’s orbit\n(Cygnus)'), (180, 'galactic anticentre\n(Auriga)'), (270, '(Vela)')]

def proj(X, Y, Z):
    # screen: right = toward l = 270°, up = toward the galactic centre (seen at a slant) and toward the north galactic pole
    return -Y, X*SE + Z*CE

def load(fname, name_col, cls_col, dist_col, extra=None):
    rows = []
    for r in csv.DictReader(open(os.path.join(here, fname), encoding='utf-8')):
        try: l, b, d = float(r['Galactic l']), float(r['Galactic b']), float(r[dist_col])
        except (ValueError, KeyError): continue
        lr, br = math.radians(l), math.radians(b)
        rows.append(dict(name=r[name_col], proper=r.get('Proper Name', r.get('Proper name', '')), cls=r[cls_col], d=d,
                         X=d*math.cos(br)*math.cos(lr), Y=d*math.cos(br)*math.sin(lr), Z=d*math.sin(br),
                         V=float(r.get('Visual Mag', r.get('V mag', 6)) or 6)))
    return rows

def plane_map(rows, R, rings, fname, title, label, note):
    fig, ax = plt.subplots(figsize=(13, 10), dpi=130); fig.patch.set_facecolor('black'); ax.set_facecolor('black')
    # the plane: rings and longitude spokes
    for rr in rings:
        ax.add_patch(Ellipse((0, 0), 2*rr, 2*rr*SE, fill=False, ec='#2f5f8f' if rr != R else '#4a86c0', lw=.9 if rr != R else 1.4))
        u, v = proj(rr*math.cos(math.radians(160)), rr*math.sin(math.radians(160)), 0)
        ax.text(u, v, '%g ly' % rr, color='#5c8dbd', fontsize=8, ha='center', va='top')
    for l in range(0, 360, 30):
        u, v = proj(R*math.cos(math.radians(l)), R*math.sin(math.radians(l)), 0)
        ax.plot([0, u], [0, v], color='#1f3f60', lw=.6, zorder=0)
        u2, v2 = proj(1.05*R*math.cos(math.radians(l)), 1.05*R*math.sin(math.radians(l)), 0)
        ax.text(u2, v2, '%d°' % l, color='#7fa8d0', fontsize=8.5, ha='center', va='center')
    for l, txt in DIRS:
        f = 1.36 if l in (90, 270) else 1.2
        u, v = proj(f*R*math.cos(math.radians(l)), f*R*math.sin(math.radians(l)), 0)
        ax.text(u, v, txt, color='#c8d8ea', fontsize=9, ha='center', va='center', style='italic')
    # stars, far side first so near stalks are drawn on top
    for r in sorted(rows, key=lambda r: -r['X']):
        u0, v0 = proj(r['X'], r['Y'], 0); u, v = proj(r['X'], r['Y'], r['Z'])
        above = r['Z'] >= 0
        ax.plot([u0, u], [v0, v], color='#8fa6bf' if above else '#6d7f93', lw=.7, ls='-' if above else (0, (2, 2)), zorder=2)
        ax.plot(u0, v0, marker='o', ms=1.6, color='#8fa6bf', zorder=2)
        c = COL.get(r['cls'][:1].upper(), '#dddddd')
        s = max(10, min(120, (7.5 - max(r['V'], -1.5))**2 * 2.2))
        ax.scatter([u], [v], s=s, color=c, edgecolors='none', zorder=3)
        lab = label(r)
        if lab: ax.annotate(lab, (u, v), xytext=(5, 3), textcoords='offset points', color='#f0f0f0', fontsize=8, zorder=4)
    us, vs = proj(0, 0, 0); ax.scatter([us], [vs], s=70, color='#ffe680', marker='*', zorder=5); ax.annotate('Sun', (us, vs), xytext=(6, -10), textcoords='offset points', color='#ffe680', fontsize=9)
    ax.set_xlim(-1.55*R, 1.55*R); ax.set_ylim(-1.32*R, 1.05*R); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, color='white', fontsize=15, pad=10)
    # legends: colours and stalks
    y0 = -1.18*R
    for k, (t, c) in enumerate([('O/B', COL['B']), ('A', COL['A']), ('F', COL['F']), ('G', COL['G']), ('K', COL['K']), ('M', COL['M']), ('white dwarf', COL['D'])]):
        ax.scatter([-1.5*R + k*0.24*R], [y0], s=40, color=c); ax.text(-1.5*R + k*0.24*R + .03*R, y0, t, color='#ccc', fontsize=8, va='center')
    ax.text(-1.5*R, y0 - .07*R, 'Star colours by spectral type.  Solid stalk: above the galactic plane.  Dashed stalk: below it.  Stalk length = height above or below the plane.',
            color='#aaa', fontsize=8)
    ax.text(-1.5*R, y0 - .11*R, note, color='#888', fontsize=7.5)
    fig.savefig(os.path.join(here, fname), bbox_inches='tight', facecolor='black'); plt.close(fig)

ten = load('stars_within_10ly.csv', 'Catalog Name', 'Classification', 'Dist (ly)')
fifty = load('stars_within_50ly.csv', 'Catalog Name', 'Classification', 'Dist (ly)')
big = load('naked_eye_stars_within_250ly.csv', 'Name', 'Spectral type', 'Distance (ly)')

TEN_LAB = {'Alpha Centauri A': 'Alpha Centauri A, B', 'Proxima Centauri': 'Proxima', 'Alpha Canis Majoris A': 'Sirius A, B',
           'Luyten 726-8 A': 'UV & BL Ceti'}
plane_map(ten, 10, [2.5, 5, 7.5, 10], 'plane_map_within_10ly.png', 'The stars within 10 light years',
          lambda r: TEN_LAB.get(r['name'], '' if r['name'] in ('Alpha Centauri B', 'Alpha Canis Majoris B', 'Luyten 726-8 B') else r['name']),
          'Hipparcos 2007 and Gaia DR3. Galactic longitude is marked around the plane.')
plane_map(fifty, 50, [10, 20, 30, 40, 50], 'plane_map_within_50ly.png', 'Naked-eye stars within 50 light years',
          lambda r: r['proper'] if r['proper'] and r['proper'] != 'Sol' else ('Zeta Reticuli' if r['name'] == 'Zeta Reticuli' else ('' if r['V'] > 3.3 else r['name'])),
          'Hipparcos 2007; layout after the Atlas of the Universe. Galactic longitude is marked around the plane.')
B250 = [r for r in big if r['d'] <= 250 and r['V'] <= 4.0]
plane_map(B250, 250, [50, 100, 150, 200, 250], 'plane_map_within_250ly.png', 'Stars within 250 light years brighter than magnitude 4 (%d)' % len(B250),
          lambda r: {'Alpha CMa': 'Sirius', 'Alpha Boo': 'Arcturus', 'Alpha Lyr': 'Vega', 'Alpha Aur': 'Capella', 'Alpha Tau': 'Aldebaran', 'Alpha Eri': 'Achernar',
                     'Alpha Vir': 'Spica', 'Alpha Leo': 'Regulus', 'Alpha PsA': 'Fomalhaut', 'Beta UMa': 'Merak', 'Zeta UMa': 'Mizar', 'Alpha Aql': 'Altair',
                     'Alpha CMi': 'Procyon', 'Beta Gem': 'Pollux', 'Alpha UMa': 'Dubhe', 'Eta UMa': 'Alkaid', 'Alpha Gem': 'Castor', 'Beta Leo': 'Denebola', 'Alpha Oph': 'Rasalhague', 'Alpha Cen': 'Alpha Centauri', 'Alpha Ari': 'Hamal', 'Alpha CrB': 'Alphecca', 'Beta Cas': 'Caph', 'Alpha Hya': 'Alphard', 'Gamma Dra': 'Eltanin'}.get(r['name'], '') if r['d'] > 60 else '',
          'Hipparcos 2007, stars of magnitude 4.0 or brighter; named stars labelled beyond 60 ly (the nearer ones are in the 50 ly map). Galactic longitude is marked around the plane.')
print('done')
