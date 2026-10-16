# Close-up (±45 ly) of the Sirius A-B orbit frame: nearby supercluster members and well-known non-member stars.
import math, os, pickle
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'; RC = '#ffe24a'
D = pickle.load(open(os.path.join(here, 'orbit_frame.pkl'), 'rb'))
mem = [s for s in D['out'] if s['grp'] == 'stream' and s['dist'] <= 47]
REFS = pickle.load(open(os.path.join(here, 'reference_stars_frame.pkl'), 'rb'))
OFF = {1: {'Sun': (8, -4), 'Alpha Centauri': (7, 2), 'Barnard’s Star': (7, 2), 'Wolf 359': (-56, 4), 'Procyon': (7, -4), 'Epsilon Eridani': (7, -4), 'Tau Ceti': (7, 2),
           '61 Cygni': (7, 2), 'Altair': (7, 2), 'Vega': (7, 2), 'Fomalhaut': (7, 2), 'Pollux': (7, 2), 'Capella': (7, 2), 'Arcturus': (7, 2)},
       2: {'Sun': (8, -4), 'Alpha Centauri': (-94, -7), 'Barnard’s Star': (7, -4), 'Wolf 359': (-58, 5), 'Procyon': (-58, 6), 'Epsilon Eridani': (7, -4), 'Tau Ceti': (7, 2),
           '61 Cygni': (-54, 4), 'Altair': (-40, -9), 'Vega': (-34, 3), 'Fomalhaut': (7, 2), 'Pollux': (-40, 4), 'Capella': (7, 2), 'Arcturus': (7, 2)}}
MOFF = {1: {}, 2: {'χ¹ Ori': (9, -2)}}
fig = plt.figure(figsize=(20, 10.5), dpi=125); fig.patch.set_facecolor(BG)
for k, (i, j, tt, xl, yl) in enumerate(((0, 1, 'Looking down on the Sirius A–B orbit plane from its pole', 'X: light years toward Sirius B’s periastron', 'Y: light years in the direction B moves'),
                                        (0, 2, 'Edge-on to the orbit plane', 'X: light years toward Sirius B’s periastron', 'Z: light years toward the orbit pole'))):
    a_ = fig.add_axes([0.04 + 0.49*k, 0.1, 0.44, 0.8]); a_.set_facecolor(BG); a_.tick_params(colors='#9aa', labelsize=10); a_.grid(color='#1c2636', lw=.6)
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
    if j == 2: a_.axhline(0, color='#ffd27a', lw=.9, alpha=.6); a_.text(-44, 0.8, 'orbit plane', color='#ffd27a', fontsize=9.5)
    for s in mem:
        a_.scatter([s['d'][i]], [s['d'][j]], s=max(14, (7.5 - s['V'])**2*7), color='#7fb2ff', zorder=4)
        a_.annotate(s['name'], (s['d'][i], s['d'][j]), xytext=MOFF[j].get(s['name'], (7, 3)), textcoords='offset points', color='#cfe0ff', fontsize=11)
    for r in REFS:
        a_.scatter([r['d'][i]], [r['d'][j]], s=200 if r['name'] == 'Sun' else 46, facecolors=RC if r['name'] == 'Sun' else 'none', edgecolors=RC, marker='*' if r['name'] == 'Sun' else 'o', lw=1.3, zorder=5)
        a_.annotate(r['name'], (r['d'][i], r['d'][j]), xytext=OFF[j].get(r['name'], (7, 2)), textcoords='offset points', color=RC, fontsize=10.5)
    a_.scatter([0], [0], s=200, color='white', zorder=6); a_.annotate('Sirius A and B', (0, 0), xytext=(8, 10), textcoords='offset points', color='white', fontsize=11.5, weight='bold')
    a_.set_xlim(-45, 45); a_.set_ylim(-45, 45); a_.set_aspect('equal')
    a_.set_xlabel(xl, color='#ccc', fontsize=11); a_.set_ylabel(yl, color='#ccc', fontsize=11); a_.set_title(tt, color='white', fontsize=13)
fig.suptitle('Around Sirius, in the frame of its A–B orbit: nearby supercluster members (blue) and well-known non-members (yellow)', color='white', fontsize=15)
fig.text(.04, .02, 'Same frame as the supercluster charts: origin at Sirius, Z along the pole of the Sirius B orbit, X toward periastron. Close-up of ±45 light years. '
         'The Big Dipper nucleus lies off this chart, about 84 ly away toward X = −72, Y = −35, Z = −28. Positions: Hipparcos and the Stars-within-10/50-ly lists.', color='#9aa', fontsize=10, wrap=True)
fig.savefig(os.path.join(here, 'sirius_neighbourhood_landmarks_close_up.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print('done')
