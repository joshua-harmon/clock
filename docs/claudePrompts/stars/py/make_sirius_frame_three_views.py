# Three views of the same stars in the Sirius A-B orbit frame: from the orbit pole, and the two 90° turns from it
# (tipped about the X axis: X across, Z up; tipped about the Y axis: Y across, Z up).
import math, os, pickle, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'; RC = '#ffe24a'
D = pickle.load(open(os.path.join(here, 'orbit_frame.pkl'), 'rb')); REFS = pickle.load(open(os.path.join(here, 'reference_stars_frame.pkl'), 'rb'))
mode = sys.argv[1] if len(sys.argv) > 1 else 'close'
BD = pickle.load(open(os.path.join(here, 'brown_dwarfs_frame.pkl'), 'rb')) if mode in ('bd', 'bdzoom') else []
BDC = {'L': '#ff9a3c', 'T': '#e05cd0', 'Y': '#9d7bff', 'M': '#ff9a3c'}
if mode in ('close', 'bd', 'bdzoom'):
    R = 17 if mode == 'bdzoom' else 45; mem = [s for s in D['out'] if s['grp'] == 'stream' and s['dist'] <= 47]; refs = REFS
    out_png = 'sirius_neighbourhood_three_views.png'; ttl = 'Around Sirius (±45 ly): supercluster members (blue) and well-known non-members (yellow), three views of the Sirius A–B orbit frame'
    if mode == 'bdzoom':
        out_png = 'sirius_brown_dwarfs_close_up_three_views.png'
        ttl = 'The nearest brown dwarfs to Sirius (±17 ly) in the Sirius A–B orbit frame: L dwarfs (orange), T dwarfs (pink), Y dwarfs (violet), with nearby stars (yellow)'
    if mode == 'bd':
        out_png = 'sirius_neighbourhood_brown_dwarfs_three_views.png'
        ttl = 'Brown dwarfs around Sirius (±45 ly) in the Sirius A–B orbit frame: L dwarfs (orange), T dwarfs (pink), Y dwarfs (violet), with stars for reference'
else:
    R = 105; mem = [s for s in D['out'] if s['grp'] in ('stream', 'core') and s['dist'] <= 100]; refs = REFS
    out_png = 'sirius_supercluster_100ly_three_views.png'; ttl = 'The Sirius supercluster within 100 ly, three views of the Sirius A–B orbit frame (red: Big Dipper nucleus; blue: stream; yellow: non-members)'
OFS = {('Wolf 359', 0): (-50, 2), ('Procyon', 1): (-50, 5), ('Pollux', 1): (-40, 4), ('Wolf 359', 1): (-52, -9), ('61 Cygni', 1): (-50, 4), ('Altair', 1): (-38, -10), ('Barnard’s Star', 1): (6, -5),
       ('Procyon', 2): (-50, 4), ('Epsilon Eridani', 2): (-82, -5), ('Tau Ceti', 2): (-50, -9), ('Sun', 2): (8, -2), ('Gliese 519', 2): (-62, -5)}
BOFS = {('WISE 0855−0714', 0): (-95, 8), ('2MASS 0939−2448', 0): (-110, -12),
        ('2MASS 0939−2448', 1): (-120, 6), ('Scholz’s star B', 1): (6, 6), ('WISE 0855−0714', 1): (-130, -2), ('Luhman 16 A, B', 1): (6, -12),
        ('ε Indi Ba, Bb', 1): (6, 6), ('WISE 0521+1025', 1): (6, -20), ('2MASS 0415−0935', 1): (6, 6), ('DENIS 0255−4700', 1): (6, 6),
        ('DENIS 0255−4700', 2): (6, 6), ('Luhman 16 A, B', 2): (6, -12), ('2MASS 0939−2448', 2): (-120, 8), ('2MASS 1114−2618', 2): (6, 8), ('DENIS 0817−6155', 2): (6, -12)}
mem = [s for s in mem if mode != 'bdzoom' or s['dist'] <= 16.6]
LABM = {s['name'] for s in mem if s['grp'] == 'stream' and s['name'] not in ('χ Cet B (EZ Cet)', 'HD 238224')}
NICK = {'χ Cet A': 'χ Cet A, B'} if mode != 'close' else {}
VIEWS = [((0, 1), 'From the orbit pole (looking down on the orbit plane)', 'X: toward Sirius B’s periastron (ly)', 'Y: direction Sirius B moves (ly)', None),
         ((0, 2), 'Turned 90° about the X axis (orbit plane edge-on, seen from −Y)', 'X: toward periastron (ly)', 'Z: toward the orbit pole (ly)', 'h'),
         ((1, 2), 'Turned 90° about the Y axis (orbit plane edge-on, seen from +X)', 'Y: direction Sirius B moves (ly)', 'Z: toward the orbit pole (ly)', 'h')]
fig = plt.figure(figsize=(24, 9.2), dpi=120); fig.patch.set_facecolor(BG)
for k, ((i, j), tt, xl, yl, line) in enumerate(VIEWS):
    a_ = fig.add_axes([0.03 + 0.33*k, 0.12, 0.29, 0.76]); a_.set_facecolor(BG); a_.tick_params(colors='#9aa', labelsize=9); a_.grid(color='#1c2636', lw=.6)
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
    if line: a_.axhline(0, color='#ffd27a', lw=.9, alpha=.6); a_.text(-R*.97, R*.015, 'orbit plane', color='#ffd27a', fontsize=9)
    else: a_.add_patch(plt.Circle((0, 0), R*.04, fill=False, ec='#ffd27a', lw=.8, alpha=.6))
    for s in sorted(mem, key=lambda s: -s['V']):
        c = '#ff6b5e' if s['grp'] == 'core' else '#7fb2ff'
        a_.scatter([s['d'][i]], [s['d'][j]], s=max(10, (7.5 - s['V'])**2*6), color=c, zorder=4)
        if s['name'] in LABM and s['name'] != 'HD 238224':
            a_.annotate(NICK.get(s['name'], s['name']), (s['d'][i], s['d'][j]), xytext=OFS.get((s['name'], k), (6, 3)), textcoords='offset points', color='#ffb0a8' if s['grp'] == 'core' else '#cfe0ff', fontsize=9)
    for r in refs:
        if max(abs(r['d'][i]), abs(r['d'][j])) > R or (mode == 'bdzoom' and r['dist'] > 16.6): continue
        a_.scatter([r['d'][i]], [r['d'][j]], s=170 if r['name'] == 'Sun' else 38, facecolors=RC if r['name'] == 'Sun' else 'none', edgecolors=RC, marker='*' if r['name'] == 'Sun' else 'o', lw=1.2, zorder=5)
        if mode in ('close', 'bd', 'bdzoom'):
            a_.annotate(r['name'], (r['d'][i], r['d'][j]), xytext=OFS.get((r['name'], k), (6, -10) if r['name'] == 'Sun' else (6, 2)), textcoords='offset points', color=RC, fontsize=9)
    if mode == 'wide':
        cc = D['cc']; a_.annotate('Big Dipper\nnucleus', (cc[i], cc[j]), xytext=(-30, -42), textcoords='offset points', color='#ff8f80', fontsize=10, ha='center',
                                   arrowprops=dict(arrowstyle='-', color='#ff8f80', lw=.7))
    for b in BD:
        if max(abs(b['d'][i]), abs(b['d'][j])) > R or (mode == 'bdzoom' and b['dist'] > 16.6): continue
        big = mode == 'bdzoom'
        a_.scatter([b['d'][i]], [b['d'][j]], s=60 if big else 16, marker='D', color=BDC[b['cls']], alpha=.95 if big else .7, zorder=3, edgecolors='none')
        if big:
            a_.annotate((b['name'] or 'L/T dwarf') + ' (' + b['spt'] + ')', (b['d'][i], b['d'][j]), xytext=BOFS.get((b['name'], k), (6, -10)), textcoords='offset points', color=BDC[b['cls']], fontsize=9)
    a_.scatter([0], [0], s=170, color='white', zorder=6); a_.annotate('Sirius A, B', (0, 0), xytext=(7, 7), textcoords='offset points', color='white', fontsize=10, weight='bold')
    a_.set_xlim(-R, R); a_.set_ylim(-R, R); a_.set_aspect('equal')
    a_.set_xlabel(xl, color='#ccc'); a_.set_ylabel(yl, color='#ccc'); a_.set_title(tt, color='white', fontsize=12)
fig.suptitle(ttl, color='white', fontsize=14.5)
foot = {'bd': 'Brown dwarfs: the 20-pc census of L, T and Y dwarfs (Kirkpatrick et al. 2021), 227 systems within 47 ly of Sirius; names of the nearest are in the close-up chart. ',
        'bdzoom': 'Shown: every brown dwarf and the listed stars within 16.6 ly of Sirius. Brown dwarfs: the 20-pc census of L, T and Y dwarfs (Kirkpatrick et al. 2021); binaries are one point with both types. '}.get(mode, '')
fig.text(.03, .02, foot + 'Frame: origin at Sirius; X toward Sirius B’s periastron, Y the direction B moves there, Z the orbit pole. Left: looking down the pole. Middle: the same scene tipped 90° about the X axis. '
         'Right: tipped 90° about the Y axis, looking back from the periastron side. The front/back sense of the orbit tilt follows the usual radial-velocity convention.' +
         (' Non-member labels are given in the close-up chart.' if mode == 'wide' else ''), color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, out_png), facecolor=BG, bbox_inches='tight'); plt.close(fig); print('done', out_png)
