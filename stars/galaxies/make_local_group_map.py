# The Local Group in three dimensions: galaxies within 1.8 Mpc from the Updated Nearby Galaxy Catalog
# (Karachentsev, Makarov & Kaisina 2013, VizieR J/AJ/145/101). Sun-centred galactic coordinates, million light years.
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from binary_stars import T
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'; MLY = 3.26156     # Mpc -> million light years
gal = []
for item in open(os.path.join(here, 'local_group_karachentsev2013.txt')).read().replace('\n', ' | ').split('|'):
    f = item.split()
    if not f: continue
    nm, ra, de, d, bm, md = f[0], float(f[1]), float(f[2]), float(f[3]), f[4], f[5]
    if nm == 'MESSIER032': d = 0.78        # catalog lists 0.49 Mpc; M32 is M31's close companion (~0.78 Mpc)
    a, dd = math.radians(ra), math.radians(de)
    p = T @ (d*MLY*np.array([math.cos(dd)*math.cos(a), math.cos(dd)*math.sin(a), math.sin(dd)]))
    gal.append(dict(name=nm, p=p, M=float(bm) if bm != 'na' else -6.0, host=md, d=d*MLY))
print(len(gal), 'galaxies')
NICE = {'Milky_Way': 'Milky Way', 'MESSIER031': 'Andromeda (M31)', 'MESSIER033': 'Triangulum (M33)', 'MESSIER032': 'M32', 'NGC0205': 'M110',
        'LMC': 'Large Magellanic Cloud', 'SMC': 'Small Magellanic Cloud', 'NGC6822': 'NGC 6822', 'IC1613': 'IC 1613', 'IC0010': 'IC 10', 'WLM': 'WLM',
        'LeoA': 'Leo A', 'SexA': 'Sextans A', 'SexB': 'Sextans B', 'NGC3109': 'NGC 3109', 'UGC04879': 'UGC 4879', 'HIZSS003': 'HIZSS 3', 'Phoenix': 'Phoenix', 'Tucana': 'Tucana', 'Pegasus': 'Pegasus dIrr',
        'NGC0185': 'NGC 185', 'NGC0147': 'NGC 147', 'Sag_dIr': 'Sagittarius dIrr', 'DDO210': 'Aquarius (DDO 210)', 'Cetus': 'Cetus', 'LeoI': 'Leo I', 'Fornax': 'Fornax'}
OUT = {'SexA', 'SexB', 'NGC3109', 'Antlia', 'HIZSS003', 'UGC04879'}
def colour(g):
    if g['name'] in OUT: return '#9aa'
    if g['name'] == 'Milky_Way' or g['host'] in ('Milky_Way', 'LMC', 'Sag_dSph'): return '#7fb2ff'
    if g['name'] in OUT: return '#9aa'
    if g['name'] == 'MESSIER033' or g['host'] == 'MESSIER033': return '#ffb36b'
    return '#ff6b5e'
def size(M): return max(6, (M + 2)**2*1.6) if M < -2 else 6
m31 = next(g for g in gal if g['name'] == 'MESSIER031')['p']; mw = next(g for g in gal if g['name'] == 'Milky_Way')['p']
bary = (mw*1.0 + m31*1.3)/2.3                      # rough barycentre, M31 taken ~1.3 x the Milky Way's mass
print('MW-M31 separation %.2f Mly; M31 l/b' % np.linalg.norm(m31 - mw), np.degrees(math.atan2(m31[1], m31[0])) % 360, np.degrees(math.asin(m31[2]/np.linalg.norm(m31))))
fig = plt.figure(figsize=(22, 10.5), dpi=120); fig.patch.set_facecolor(BG)
VIEWS = [((0, 1), 'Seen from above the Milky Way’s disc (north galactic pole)', 'million ly toward the galactic centre →', 'million ly toward l = 90° (Cygnus)'),
         ((1, 2), 'Turned 90°: seen from the galactic centre direction', 'million ly toward l = 90° →', 'million ly north of the galactic plane')]
OFF = {'Small Magellanic Cloud': (8, -12), 'M32': (8, -12), 'M110': (8, 6), 'NGC 147': (-50, 6), 'NGC 185': (8, 6)}
for k, ((i, j), tt, xl, yl) in enumerate(VIEWS):
    ax = fig.add_axes([0.04 + 0.49*k, 0.1, 0.44, 0.8]); ax.set_facecolor(BG); ax.tick_params(colors='#9aa'); ax.grid(color='#1c2636', lw=.6)
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.add_patch(plt.Circle((mw[i], mw[j]), 0.05, fill=False, ec='#7fb2ff', lw=.8, alpha=.7))
    ax.add_patch(plt.Circle((m31[i], m31[j]), 0.11, fill=False, ec='#ff6b5e', lw=.8, alpha=.7))
    for R in (1.0, 2.0, 3.0): ax.add_patch(plt.Circle((bary[i], bary[j]), R, fill=False, ec='#3a4a60', lw=.7, ls=(0, (3, 4))))
    ax.text(bary[i] + 2.05*0.71, bary[j] + 2.05*0.71, '2 Mly from the\ncentre of mass', color='#667', fontsize=8.5)
    ax.scatter([bary[i]], [bary[j]], marker='+', s=90, color='#ffe24a', zorder=5)
    for g in sorted(gal, key=lambda g: -g['M']):
        ax.scatter([g['p'][i]], [g['p'][j]], s=size(g['M']), color=colour(g), alpha=.9, edgecolors='none', zorder=4)
        nm = NICE.get(g['name'])
        if g['name'] == 'NGC3109': nm = 'NGC 3109 + Antlia'
        if nm and (g['name'] in OUT or g['M'] < -11 or nm in ('Leo I', 'Fornax')) and not (k == 0 and nm in ('Leo I', 'Fornax')):
            ax.annotate(nm, (g['p'][i], g['p'][j]), xytext=OFF.get(nm, (7, 4)), textcoords='offset points', color='#e6e6e6', fontsize=9.5 if g['M'] < -17 else 8.5,
                        weight='bold' if g['M'] < -18 else 'normal')
    ax.set_xlim(-4.8, 2.6); ax.set_ylim(-4.4, 4.4); ax.set_aspect('equal')
    if k == 1: ax.set_xlim(-4.6, 3.8); ax.set_ylim(-3.2, 4.0)
    ax.set_xlabel(xl, color='#ccc'); ax.set_ylabel(yl, color='#ccc'); ax.set_title(tt, color='white', fontsize=13)
fig.suptitle('The Local Group: the Milky Way, Andromeda and their companions (82 galaxies within 5.9 million light years)', color='white', fontsize=16)
fig.text(.04, .02, 'Blue: the Milky Way and its satellites. Red: Andromeda (M31) and its satellites. Orange: Triangulum (M33). Grey: outlying galaxies just beyond the Local Group proper (the NGC 3109 – Antlia – Sextans A/B group, UGC 4879 and HIZSS 3). '
         'Dot size by brightness; small circles show the sizes of the two big discs. Yellow +: rough centre of mass. Distances: Updated Nearby Galaxy Catalog (Karachentsev et al. 2013); '
         'M32 placed at Andromeda’s distance.', color='#9aa', fontsize=10, wrap=True)
fig.savefig(os.path.join(here, 'local_group_map.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig); print('done')
