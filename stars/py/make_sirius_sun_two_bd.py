# Sirius, the Sun and the two nearest brown dwarfs (Luhman 16 A,B nearest the Sun; WISE 0855-0714 nearest Sirius)
# in the Sirius A-B orbit frame, three views. Dotted lines are labelled with true 3-D distances.
import math, os, pickle, itertools
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
BD = {b['name']: b for b in pickle.load(open(os.path.join(here, 'brown_dwarfs_frame.pkl'), 'rb'))}
REF = {r['name']: r for r in pickle.load(open(os.path.join(here, 'reference_stars_frame.pkl'), 'rb'))}
P = [('Sirius A, B', np.zeros(3), 'white', 'o', 260, 'A1V + white dwarf'),
     ('Sun', REF['Sun']['d'], '#ffe24a', '*', 380, 'G2V'),
     ('Luhman 16 A, B', BD['Luhman 16 A, B']['d'], '#ff9a3c', 'D', 120, 'L7.5 + T0.5 brown dwarfs'),
     ('WISE 0855−0714', BD['WISE 0855−0714']['d'], '#9d7bff', 'D', 120, 'Y4 brown dwarf')]
for a, b in itertools.combinations(P, 2): print('%s - %s: %.1f ly' % (a[0], b[0], np.linalg.norm(a[1] - b[1])))
VIEWS = [((0, 1), 'From the orbit pole (looking down on the orbit plane)', 'X: toward Sirius B’s periastron (ly)', 'Y: direction Sirius B moves (ly)', False),
         ((0, 2), 'Turned 90° about the X axis (orbit plane edge-on)', 'X: toward periastron (ly)', 'Z: toward the orbit pole (ly)', True),
         ((1, 2), 'Turned 90° about the Y axis (orbit plane edge-on)', 'Y: direction Sirius B moves (ly)', 'Z: toward the orbit pole (ly)', True)]
LOFF = {('Luhman 16 A, B', 0): (10, 8), ('WISE 0855−0714', 0): (-12, -18), ('Sun', 0): (12, 4), ('Sirius A, B', 0): (10, -16),
        ('Luhman 16 A, B', 1): (-6, 30), ('WISE 0855−0714', 1): (-12, -26), ('Sun', 1): (12, -4), ('Sirius A, B', 1): (12, 22),
        ('Luhman 16 A, B', 2): (-8, 30), ('WISE 0855−0714', 2): (12, -26), ('Sun', 2): (12, -4), ('Sirius A, B', 2): (-14, 30)}
fig = plt.figure(figsize=(24, 9), dpi=120); fig.patch.set_facecolor(BG)
R = 10
for k, ((i, j), tt, xl, yl, line) in enumerate(VIEWS):
    a_ = fig.add_axes([0.03 + 0.33*k, 0.13, 0.29, 0.75]); a_.set_facecolor(BG); a_.tick_params(colors='#9aa', labelsize=10); a_.grid(color='#1c2636', lw=.6)
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
    if line: a_.axhline(0, color='#ffd27a', lw=.9, alpha=.6); a_.text(-R*.96, .15, 'orbit plane', color='#ffd27a', fontsize=9.5)
    for a, b in itertools.combinations(P, 2):
        a_.plot([a[1][i], b[1][i]], [a[1][j], b[1][j]], color='#556079', lw=.9, ls=':', zorder=1)
        if k: continue                       # distances are labelled in the open, top view only
        f = {('Sun', 'WISE 0855−0714'): .62, ('Luhman 16 A, B', 'WISE 0855−0714'): .35}.get((a[0], b[0]), .5)
        mx, my = a[1][i] + f*(b[1][i] - a[1][i]), a[1][j] + f*(b[1][j] - a[1][j])
        a_.text(mx, my, '%.1f ly' % np.linalg.norm(a[1] - b[1]), color='#8fa0bd', fontsize=9, ha='center', va='center', zorder=2,
                bbox=dict(boxstyle='round,pad=0.15', fc=BG, ec='none'))
    for nm, v, c, mk, sz, sp_ in P:
        a_.scatter([v[i]], [v[j]], s=sz, color=c, marker=mk, zorder=5, edgecolors='none')
        o = LOFF.get((nm, k), (10, 6))
        a_.annotate('%s\n%s' % (nm, sp_), (v[i], v[j]), xytext=o, textcoords='offset points', color=c, fontsize=10.5, weight='bold' if nm == 'Sirius A, B' else 'normal',
                    ha='left' if o[0] > 0 else 'right', va='center')
    a_.set_xlim(-R, R); a_.set_ylim(-R, R); a_.set_aspect('equal')
    a_.set_xlabel(xl, color='#ccc', fontsize=11); a_.set_ylabel(yl, color='#ccc', fontsize=11); a_.set_title(tt, color='white', fontsize=12.5)
fig.suptitle('Sirius, the Sun and the two nearest brown dwarfs: Luhman 16 (nearest the Sun, 6.5 ly) and WISE 0855−0714 (nearest Sirius, 4.7 ly)', color='white', fontsize=15)
fig.text(.03, .025, 'Frame: origin at Sirius; X toward Sirius B’s periastron, Y the direction B moves there, Z the orbit pole. Dotted lines join each pair and are labelled with the true '
         '3-D distance, written on the left view. Brown dwarfs: Kirkpatrick et al. (2021) 20-pc census. The front/back sense of the orbit tilt follows the usual radial-velocity convention.',
         color='#9aa', fontsize=10, wrap=True)
fig.savefig(os.path.join(here, 'sirius_sun_and_nearest_brown_dwarfs.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig); print('done')
