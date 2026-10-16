# Charts for make_sirius_frame_uma.py: the Ursa Major group seen in the frame of the Sirius A-B orbit.
import math, os, pickle
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
D = pickle.load(open(os.path.join(here, 'orbit_frame.pkl'), 'rb')); out, M = D['out'], D['M']
sun_lam, sun_bet, sunv = D['sun']; cc = D['cc']
COL = {'core': '#ff6b5e', 'stream': '#7fb2ff', 'sirius': '#ffffff'}
def sz(V): return max(6, (7.5 - V)**2*5) if V < 7.5 else 6
def wrap(l): return -math.radians(((l + 180) % 360) - 180)       # longitude increasing to the left, 0° in the centre
LABEL = {'Menkalinan (β Aur)', 'Alphecca (α CrB)', 'Skat (δ Aqr)', 'Adhafera (ζ Leo)', 'γ Lep', 'β Ser', 'ζ Boo', 'χ¹ Ori', '21 LMi', 'χ Cet A', 'γ Mic', 'ζ Crt',
         'ζ TrA', '16 Lyr', '66 Tau', '59 Dra', '89 Psc', '18 Boo', 'π¹ UMa', '29 Com', '6 Sex', '41 Vir'}
fig = plt.figure(figsize=(21, 18), dpi=120); fig.patch.set_facecolor(BG)
ax = fig.add_axes([0.03, 0.55, 0.94, 0.38], projection='hammer'); ax.set_facecolor(BG)
ax.grid(color='#2a3548', lw=.6); ax.tick_params(colors='#9aa', labelsize=9)
ax.set_xticks([math.radians(x) for x in range(-150, 181, 30)]); ax.set_xticklabels(['%d°' % ((-x) % 360) for x in range(-150, 181, 30)], color='#9aa')
# the orbit plane (equator), the galactic plane and the galactic poles
ax.plot(np.linspace(-math.pi, math.pi, 200), np.zeros(200), color='#ffd27a', lw=1.6)
ax.text(wrap(25), math.radians(2.5), 'plane of the Sirius A–B orbit', color='#ffd27a', fontsize=10)
g = np.radians(np.arange(0, 360.5, 0.5)); gp = np.array([M @ np.array([math.cos(t), math.sin(t), 0]) for t in g])
gl = np.array([wrap(math.degrees(math.atan2(v[1], v[0]))) for v in gp]); gb = np.arcsin(gp[:, 2])
brk = np.where(np.abs(np.diff(gl)) > 3)[0] + 1
for seg_l, seg_b in zip(np.split(gl, brk), np.split(gb, brk)): ax.plot(seg_l, seg_b, color='#6b8f6b', lw=1.2, ls=(0, (5, 3)))
k = np.argmin(np.abs(gb - math.radians(40))); ax.text(gl[k] + .05, gb[k], 'galactic plane', color='#8fbf8f', fontsize=10)
for v, nm in ((M @ np.array([0, 0, 1.0]), 'north galactic pole'), (M @ np.array([0, 0, -1.0]), 'south galactic pole'), (M @ np.array([1.0, 0, 0]), 'toward the galactic centre')):
    L, B = wrap(math.degrees(math.atan2(v[1], v[0]))), math.asin(v[2]); ax.scatter([L], [B], marker='x', s=50, color='#8fbf8f', zorder=4)
    ax.annotate(nm, (L, B), xytext=(6, 5), textcoords='offset points', color='#8fbf8f', fontsize=9)
ax.scatter([wrap(sun_lam)], [math.radians(sun_bet)], s=160, color='#ffe24a', marker='*', zorder=6)
ax.annotate('the Sun (8.6 ly)', (wrap(sun_lam), math.radians(sun_bet)), xytext=(8, -12), textcoords='offset points', color='#ffe24a', fontsize=10.5, weight='bold')
for s in sorted(out, key=lambda s: -s['V']):
    if s['grp'] == 'sirius': continue
    L, B = wrap(s['lam']), math.radians(s['bet'])
    ax.scatter([L], [B], s=sz(s['V']), color=COL[s['grp']], edgecolors='none', alpha=.95, zorder=5)
    if s['name'] in LABEL: ax.annotate(s['name'], (L, B), xytext=(5, 4), textcoords='offset points', color='#cfd8e6', fontsize=8.5)
L, B = wrap(math.degrees(math.atan2(cc[1], cc[0]))), math.asin(cc[2]/np.linalg.norm(cc))
ax.annotate('Big Dipper nucleus\n(12 core stars, 84 ly from Sirius)', (L, B), xytext=(30, 30), textcoords='offset points', color='#ff8f80', fontsize=10,
            arrowprops=dict(arrowstyle='-', color='#ff8f80', lw=.8))
ax.set_title('The sky seen from Sirius, in coordinates fixed to the Sirius A–B orbit (orbit plane = equator; longitude from Sirius B’s periastron, in the direction B moves)',
             color='white', fontsize=12.5, pad=14)

# lower panels: looking down the orbit pole, and edge-on
bx = fig.add_axes([0.04, 0.07, 0.44, 0.42]); cx = fig.add_axes([0.53, 0.07, 0.44, 0.42])
for a_, (i, j), tt, xl, yl in ((bx, (0, 1), 'Looking down on the orbit plane from its pole (Z up out of the page)', 'X: light years toward Sirius B’s periastron', 'Y: light years in the direction B moves'),
                               (cx, (0, 2), 'Edge-on to the orbit plane', 'X: light years toward Sirius B’s periastron', 'Z: light years toward the orbit pole')):
    a_.set_facecolor(BG); a_.tick_params(colors='#9aa', labelsize=9); a_.grid(color='#1c2636', lw=.6)
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
    a_.axhline(0, color='#ffd27a', lw=.8, alpha=.5) if j == 2 else None
    for s in sorted(out, key=lambda s: -s['V']):
        a_.scatter([s['d'][i]], [s['d'][j]], s=sz(s['V']) if s['grp'] != 'sirius' else 140, color=COL[s['grp']], edgecolors='none', alpha=.95, zorder=4)
        if s['name'] in LABEL and abs(s['d'][i]) < 230 and abs(s['d'][j]) < 230: a_.annotate(s['name'], (s['d'][i], s['d'][j]), xytext=(4, 3), textcoords='offset points', color='#cfd8e6', fontsize=8)
    a_.annotate('Sirius A and B\n(orbit 0.0005 ly across)', (0, 0), xytext=(60, -62), textcoords='offset points', color='white', fontsize=9.5, ha='left', arrowprops=dict(arrowstyle='-', color='#bbb', lw=.6))
    a_.scatter([sunv[i]], [sunv[j]], s=140, marker='*', color='#ffe24a', zorder=6); a_.annotate('Sun', (sunv[i], sunv[j]), xytext=(6, 4), textcoords='offset points', color='#ffe24a', fontsize=10)
    a_.annotate('nucleus', (cc[i], cc[j]), xytext=(-10, 16), textcoords='offset points', color='#ff8f80', fontsize=10, ha='right')
    a_.set_xlim(-240, 240); a_.set_ylim(-240, 240); a_.set_aspect('equal')
    a_.set_xlabel(xl, color='#ccc'); a_.set_ylabel(yl, color='#ccc'); a_.set_title(tt, color='white', fontsize=12)
bx.text(-232, -232, 'red: core (Big Dipper nucleus)   blue: stream members\nfar members (66 Tau, 29 Com, HD 27022, HD 220096) lie beyond 240 ly', color='#9aa', fontsize=9, va='bottom')
fig.suptitle('The Ursa Major moving group (“Sirius supercluster”) in the frame of the Sirius A–B orbit', color='white', fontsize=16)
fig.text(.03, .012, 'Frame: origin at Sirius; Z along the pole of the Sirius B orbit (Bond et al. 2017 elements; pole at galactic l = 268°, b = +6°, so the orbit plane stands 84° to the galactic plane); X toward periastron. '
         'The front/back sense of the orbit tilt follows the usual radial-velocity convention; if reversed, the frame is mirrored. Members: core after King et al. (2003) and stream stars as listed by Wikipedia; '
         'positions from Hipparcos (2007). Sirius’s own membership is disputed (King et al. find it younger than the group).', color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, 'uma_group_in_sirius_orbit_frame.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print('done')
