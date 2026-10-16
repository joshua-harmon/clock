# Charts for make_sirius_frame_uma.py: the Ursa Major group seen in the frame of the Sirius A-B orbit.
import math, os, pickle
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
D = pickle.load(open(os.path.join(here, 'orbit_frame.pkl'), 'rb')); M = D['M']
REFS = pickle.load(open(os.path.join(here, 'reference_stars_frame.pkl'), 'rb'))
RC = '#ffe24a'
out = [s for s in D['out'] if s['grp'] in ('core', 'sirius') or s['dist'] <= 100]   # the Big Dipper nucleus plus members within 100 ly of Sirius
sun_lam, sun_bet, sunv = D['sun']; cc = D['cc']
COL = {'core': '#ff6b5e', 'stream': '#7fb2ff', 'sirius': '#ffffff'}
def sz(V): return max(6, (7.5 - V)**2*5) if V < 7.5 else 6
def wrap(l): return -math.radians(((l + 180) % 360) - 180)       # longitude increasing to the left, 0° in the centre
LABEL0 = {'Menkalinan (β Aur)', 'Alphecca (α CrB)', 'Skat (δ Aqr)', 'Adhafera (ζ Leo)', 'γ Lep', 'β Ser', 'ζ Boo', 'χ¹ Ori', '21 LMi', 'χ Cet A', 'γ Mic', 'ζ Crt',
         'ζ TrA', '16 Lyr', '66 Tau', '59 Dra', '89 Psc', '18 Boo', 'π¹ UMa', '29 Com', '6 Sex', '41 Vir'}
LABEL = {s['name'] for s in out if s['grp'] == 'stream'}
CORELAB = {'Alioth', 'Mizar', 'Merak', 'Phecda', 'Megrez', 'Alcor', '37 UMa', '78 UMa'}
HOFF = {'Barnard’s Star': (6, -12), 'Sun': (-26, 5), 'Alpha Centauri': (6, 3), 'Wolf 359': (6, 3), 'Arcturus': (6, -12), 'Epsilon Eridani': (6, 3), 'Tau Ceti': (6, -10), 'Pollux': (-34, 6), 'Procyon': (6, 3)}
RPOFF = {('Sun', 1): (6, 3), ('Alpha Centauri', 1): (6, 3), ('Barnard’s Star', 1): (6, 3), ('Wolf 359', 1): (-50, 3), ('Procyon', 1): (6, -10), ('Epsilon Eridani', 1): (6, 3), ('Tau Ceti', 1): (6, 3), ('61 Cygni', 1): (-48, 3), ('Altair', 1): (6, 3), ('Vega', 1): (-34, 3),
         ('Sun', 2): (6, -10), ('Alpha Centauri', 2): (6, 3), ('Barnard’s Star', 2): (6, -10), ('Wolf 359', 2): (-52, -3), ('Procyon', 2): (6, 3), ('Epsilon Eridani', 2): (6, -10), ('Tau Ceti', 2): (6, -4), ('61 Cygni', 2): (-48, 3), ('Altair', 2): (6, -6), ('Vega', 2): (-34, 3)}
OFF = {('78 UMa', 1): (6, 6), ('Alioth', 1): (-40, -14), ('Mizar', 1): (-38, 2), ('Merak', 2): (6, 8), ('Phecda', 2): (-50, 6), ('Mizar', 2): (-40, -4), ('Megrez', 1): (8, 2)}
SOFF = {('χ Cet B (EZ Cet)', 1): (6, -12), ('χ Cet B (EZ Cet)', 2): (6, -12), ('HD 165185', 2): (8, -16)}
EDGE = {'Merak', 'Phecda', 'Mizar'}
fig = plt.figure(figsize=(21, 18), dpi=120); fig.patch.set_facecolor(BG)
ax = fig.add_axes([0.03, 0.55, 0.94, 0.38], projection='hammer'); ax.set_facecolor(BG)
ax.grid(color='#2a3548', lw=.6); ax.tick_params(colors='#9aa', labelsize=9)
ax.set_xticks([math.radians(x) for x in range(-150, 181, 30)]); ax.set_xticklabels(['%d°' % ((-x) % 360) for x in range(-150, 181, 30)], color='#9aa')
# the orbit plane (equator), the galactic plane and the galactic poles
ax.plot(np.linspace(-math.pi, math.pi, 200), np.zeros(200), color='#ffd27a', lw=1.6)
ax.text(wrap(25), math.radians(2.5), 'plane of the Sirius A–B orbit', color='#ffd27a', fontsize=10)
for s in sorted(out, key=lambda s: -s['V']):
    if s['grp'] == 'sirius': continue
    L, B = wrap(s['lam']), math.radians(s['bet'])
    ax.scatter([L], [B], s=sz(s['V']), color=COL[s['grp']], edgecolors='none', alpha=.95, zorder=5)
    if s['name'] in LABEL and s['name'] not in ('χ Cet B (EZ Cet)', 'HD 238224'):
        ax.annotate({'χ Cet A': 'χ Cet A and B'}.get(s['name'], s['name']), (L, B), xytext={'DX Lyn': (5, -11)}.get(s['name'], (5, 4)), textcoords='offset points', color='#cfd8e6', fontsize=8.5)
for s in REFS:
    L, B = wrap(s['lam']), math.radians(s['bet'])
    ax.scatter([L], [B], s=70 if s['name'] == 'Sun' else 34, facecolors='none' if s['name'] != 'Sun' else RC, edgecolors=RC, marker='*' if s['name'] == 'Sun' else 'o', lw=1.2, zorder=6)
    ax.annotate(s['name'], (L, B), xytext=HOFF.get(s['name'], (6, -9)), textcoords='offset points', color=RC, fontsize=8.5)
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
        if s['name'] in LABEL and s['name'] != 'HD 238224': a_.annotate(s['name'], (s['d'][i], s['d'][j]), xytext=SOFF.get((s['name'], j), (6, 4)), textcoords='offset points', color='#cfd8e6', fontsize=9.5)
        elif s['name'] in CORELAB and (j == 1 or s['name'] in EDGE): a_.annotate(s['name'], (s['d'][i], s['d'][j]), xytext=OFF.get((s['name'], j), (6, 3)), textcoords='offset points', color='#ffb0a8', fontsize=8.5)
    for r_ in REFS:
        a_.scatter([r_['d'][i]], [r_['d'][j]], s=110 if r_['name'] == 'Sun' else 34, facecolors=RC if r_['name'] == 'Sun' else 'none', edgecolors=RC, marker='*' if r_['name'] == 'Sun' else 'o', lw=1.1, zorder=6)
    a_.annotate('Sirius A and B\n(orbit 0.0005 ly across)', (0, 0), xytext=(70, 130), textcoords='offset points', color='white', fontsize=9.5, ha='left', arrowprops=dict(arrowstyle='-', color='#bbb', lw=.6))
    a_.set_xlim(-105, 105); a_.set_ylim(-105, 105); a_.set_aspect('equal')
    a_.set_xlabel(xl, color='#ccc'); a_.set_ylabel(yl, color='#ccc'); a_.set_title(tt, color='white', fontsize=12)
bx.text(-101, -101, 'red: Big Dipper nucleus (12 core stars)   blue: members within 100 ly of Sirius   yellow: well-known non-members (named in the close-up chart)', color='#9aa', fontsize=9.5, va='bottom')
fig.suptitle('The Sirius supercluster within 100 light years in the frame of the Sirius A–B orbit', color='white', fontsize=16)
fig.text(.03, .012, 'Shown: the Big Dipper nucleus and the supercluster (Ursa Major group) members within 100 light years of Sirius; yellow open circles mark well-known stars that are not members (the Sun, Alpha Centauri, Barnard’s Star, Procyon, Vega, ...), for reference. Frame: origin at Sirius; Z along the pole of the Sirius B orbit (Bond et al. 2017 elements); X toward periastron. '
         'The front/back sense of the orbit tilt follows the usual radial-velocity convention; if reversed, the frame is mirrored. Members: core after King et al. (2003) and stream stars as listed by Wikipedia; '
         'positions from Hipparcos (2007). Sirius’s own membership is disputed (King et al. find it younger than the group).', color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, 'sirius_supercluster_100ly_with_landmarks.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print('done')
