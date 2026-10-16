# Space velocities (U, V, W) of the naked-eye stars within 250 ly and the Ursa Major stream.
# Inputs: naked_eye_stars_within_250ly.csv (made by make_star_maps.py) and the Pulkovo compilation of radial velocities
# for Hipparcos stars (Gontcharov 2006, VizieR III/252), saved as pulkovo_radial_velocities_within_250ly.csv.
# U toward the galactic centre, V in the direction of galactic rotation, W toward the north galactic pole; km/s relative to the Sun.
# Usage: python make_stream_map.py path\to\pulkovo_radial_velocities_within_250ly.csv
import csv, math, sys, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.dirname(os.path.abspath(__file__))
rvf = sys.argv[1] if len(sys.argv) > 1 else 'pulkovo_radial_velocities_within_250ly.csv'
T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
              [ 0.4941094279, -0.4448296300,  0.7469822445],
              [-0.8676661490, -0.1980763734,  0.4559837762]])
K = 4.740470446   # km/s for 1 mas/yr of proper motion at 1 mas of parallax

RV = {}
for r in csv.DictReader(open(rvf, encoding='utf-8')):
    try: RV[int(r['HIP'])] = float(r['RV'])
    except ValueError: pass

stars = list(csv.DictReader(open(os.path.join(out, 'naked_eye_stars_within_250ly.csv'), encoding='utf-8')))
def uvw(l, b, plx, pmra, pmde, rv):
    lr, br = math.radians(l), math.radians(b)
    re = T.T @ np.array([math.cos(br)*math.cos(lr), math.cos(br)*math.sin(lr), math.sin(br)])
    ra = math.atan2(re[1], re[0]); dec = math.asin(re[2])
    p = np.array([-math.sin(ra), math.cos(ra), 0.0])
    q = np.array([-math.sin(dec)*math.cos(ra), -math.sin(dec)*math.sin(ra), math.cos(dec)])
    return T @ ((K*pmra/plx)*p + (K*pmde/plx)*q + rv*re)

res = []
for s in stars:
    h = int(s['HIP'])
    if h not in RV: continue
    U, V, W = uvw(float(s['Galactic l']), float(s['Galactic b']), float(s['Parallax (mas)']), float(s['pm RA (mas/yr)']), float(s['pm Dec (mas/yr)']), RV[h])
    res.append(dict(s=s, HIP=h, U=U, V=V, W=W, rv=RV[h], l=float(s['Galactic l']), b=float(s['Galactic b']), ly=float(s['Distance (ly)']),
                    Vmag=float(s['V mag']), core=s['Ursa Major group nucleus'] == 'yes', name=s['Proper name'] or s['Name']))

core = [r for r in res if r['core']]
U0, V0, W0 = (float(np.median([r[k] for r in core])) for k in 'UVW')
print('nucleus stars with RV:', len(core), ' median UVW = %.1f %.1f %.1f km/s' % (U0, V0, W0))
for r in core: print('  %-12s %6.1f %6.1f %6.1f' % (r['name'], r['U'], r['V'], r['W']))
for r in res: r['d'] = math.sqrt((r['U']-U0)**2 + (r['V']-V0)**2 + (r['W']-W0)**2)
CLOSE, NEAR = 3.0, 6.0
cand = [r for r in res if r['d'] <= CLOSE and not r['core']]; poss = [r for r in res if CLOSE < r['d'] <= NEAR and not r['core']]
sir = next((r for r in res if r['HIP'] == 32349), None)
print('stars with RV', len(res), ' within 3 km/s:', len(cand), ' 3-6 km/s:', len(poss))
if sir: print('Sirius UVW %.1f %.1f %.1f, differs from the nucleus by %.1f km/s' % (sir['U'], sir['V'], sir['W'], sir['d']))
print('close:', ', '.join('%s (%.0f ly)' % (r['name'], r['ly']) for r in sorted(cand, key=lambda r: r['Vmag'])[:25]))

with open(os.path.join(out, 'ursa_major_stream_velocities.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['HIP', 'Name', 'Proper name', 'Distance (ly)', 'V mag', 'Spectral type', 'Galactic l', 'Galactic b', 'Radial velocity (km/s)',
                'U (km/s)', 'V (km/s)', 'W (km/s)', 'Difference from UMa nucleus (km/s)', 'Stream match'])
    for r in sorted(res, key=lambda r: r['d']):
        s = r['s']
        w.writerow([r['HIP'], s['Name'], s['Proper name'], s['Distance (ly)'], s['V mag'], s['Spectral type'], s['Galactic l'], s['Galactic b'], r['rv'],
                    '%.1f' % r['U'], '%.1f' % r['V'], '%.1f' % r['W'], '%.1f' % r['d'],
                    'nucleus' if r['core'] else ('close (within 3 km/s)' if r['d'] <= CLOSE else ('possible (3-6 km/s)' if r['d'] <= NEAR else ''))])

fig = plt.figure(figsize=(16, 7.8), dpi=130)
ax = fig.add_subplot(1, 2, 1)
ax.scatter([r['U'] for r in res], [r['V'] for r in res], s=5, c='#9aa8b8', edgecolors='none', label='naked-eye stars within 250 ly with radial velocities (%d)' % len(res))
ax.scatter([r['U'] for r in poss], [r['V'] for r in poss], s=16, c='#f0a040', edgecolors='none', label='possible stream members, 3–6 km/s (%d)' % len(poss))
ax.scatter([r['U'] for r in cand], [r['V'] for r in cand], s=22, c='#d0302a', edgecolors='none', label='close match, within 3 km/s (%d)' % len(cand))
ax.scatter([r['U'] for r in core], [r['V'] for r in core], s=90, facecolors='none', edgecolors='k', lw=1.2, label='Big Dipper nucleus (%d)' % len(core))
if sir: ax.annotate('Sirius', (sir['U'], sir['V']), xytext=(18, -22), textcoords='offset points', fontsize=9, arrowprops=dict(arrowstyle='-', lw=.6))
ax.annotate('Ursa Major stream', (U0, V0), xytext=(U0+6, V0+14), fontsize=9.5, color='#a01c18', weight='bold')
for lab, (u, v) in {'Hyades stream': (-40, -17), 'Pleiades / Local': (-12, -22), 'Hercules stream': (-30, -50)}.items():
    ax.annotate(lab, (u, v), xytext=(u-44, v-9), fontsize=9, color='#333')
ax.set_xlim(-80, 60); ax.set_ylim(-80, 40); ax.set_aspect('equal'); ax.grid(color='#eee')
ax.set_xlabel('U: toward the galactic centre (km/s)'); ax.set_ylabel('V: direction of galactic rotation (km/s)')
ax.set_title('Velocities relative to the Sun (the U–V plane)'); ax.legend(loc='lower right', fontsize=7.5, frameon=False)
ax2 = fig.add_subplot(1, 2, 2, projection='aitoff')
X = lambda l: math.radians(-(((l + 180) % 360) - 180))
ax2.scatter([X(r['l']) for r in res], [math.radians(r['b']) for r in res], s=3, c='#c4ccd6', edgecolors='none')
ax2.scatter([X(r['l']) for r in poss], [math.radians(r['b']) for r in poss], s=[max(8, (6.6-r['Vmag'])**2*4) for r in poss], c='#f0a040', edgecolors='none')
ax2.scatter([X(r['l']) for r in cand], [math.radians(r['b']) for r in cand], s=[max(10, (6.6-r['Vmag'])**2*5) for r in cand], c='#d0302a', edgecolors='none')
ax2.scatter([X(r['l']) for r in core], [math.radians(r['b']) for r in core], s=110, facecolors='none', edgecolors='k', lw=1.2)
ax2.annotate('Big Dipper nucleus', (X(core[0]['l']), math.radians(core[0]['b'])), xytext=(10, 10), textcoords='offset points', fontsize=8)
for r in cand + poss + ([sir] if sir else []):
    if (r['s']['Proper name'] or r['Vmag'] < 2.3) and not r['core']:
        ax2.annotate(r['name'], (X(r['l']), math.radians(r['b'])), xytext=(5, 4), textcoords='offset points', fontsize=8)
ax2.set_xticklabels(['150°','120°','90°','60°','30°','0°','330°','300°','270°','240°','210°'], fontsize=8); ax2.grid(True, color='#bbb', lw=.5)
ax2.set_title('Where the stream stars are on the sky (galactic coordinates)', pad=14)
fig.suptitle('The Ursa Major stream (Eggen’s Sirius supercluster): stars moving with the Big Dipper nucleus, U, V, W = %.1f, %.1f, %.1f km/s' % (U0, V0, W0), fontsize=12.5)
fig.text(.5, .02, 'Hipparcos parallaxes and proper motions (van Leeuwen 2007); radial velocities from the Pulkovo compilation (Gontcharov 2006). Red: full 3-D velocity within 3 km/s of the nucleus; orange: 3–6 km/s.\n'
         'A matching velocity makes a star a candidate, not a proven member: some field stars share the motion by chance.', ha='center', fontsize=8, color='#555')
fig.savefig(os.path.join(out, 'ursa_major_stream_velocities.png'), bbox_inches='tight'); plt.close(fig)
print('done')
