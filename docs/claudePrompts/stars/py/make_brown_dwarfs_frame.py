# Brown dwarfs (L, T and Y dwarfs) within 47 ly of Sirius, from the 20-pc census of Kirkpatrick et al. (2021, ApJS 253, 7;
# VizieR J/ApJS/253/7/table11), placed in the Sirius A-B orbit frame. Close binaries are merged into one point.
import csv, glob, math, os, pickle
import numpy as np
from binary_stars import T
here = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(here, 'orbit_frame.pkl'), 'rb')); M = D['M']; C = -np.linalg.solve(M, D['sun'][2])
raw = []
for f in sorted(glob.glob(os.path.join(here, 'bd', 'c*.txt'))):
    for item in open(f).read().strip().split('|'):
        ra, de, plx, spt = item.split(); raw.append((float(ra), float(de), float(plx), spt))
merged = {}
for ra, de, plx, spt in raw:
    k = (round(ra, 1), round(de, 1))
    if k in merged: merged[k]['spt'] += '+' + spt
    else: merged[k] = dict(ra=ra, de=de, plx=plx, spt=spt)
NAMES = {(133.8, -7.2): 'WISE 0855−0714', (110.6, -5.7): 'UGPS 0722−0540', (162.3, -53.3): 'Luhman 16 A, B', (92.6, -21.9): 'Gliese 229 B',
         (144.9, -24.8): '2MASS 0939−2448', (63.8, -9.6): '2MASS 0415−0935', (124.4, -61.9): 'DENIS 0817−6155', (43.8, -47.0): 'DENIS 0255−4700',
         (110.0, -8.8): 'Scholz’s star B', (57.5, -57.0): 'WISE 0350−5658', (331.0, -56.8): 'ε Indi Ba, Bb', (80.4, 10.4): 'WISE 0521+1025',
         (281.3, -64.0): 'SCR 1845−6357 B', (126.3, 28.1): 'WISE 0825+2805', (62.6, 15.1): 'WISE 0410+1502', (249.9, -68.8): 'WISE 1639−6847', (168.7, -26.3): '2MASS 1114−2618', (128.9, -8.3): '2MASS 0835−0819'}
rows = []
for k, s in merged.items():
    a, d = math.radians(s['ra']), math.radians(s['de']); dsun = 3261.56/s['plx']
    p = T @ (dsun*np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])); v = M @ (p - C); dist = np.linalg.norm(v)
    cls = s['spt'][0]
    rows.append(dict(name=NAMES.get((round(s['ra'], 1), round(s['de'], 1)), ''), spt=s['spt'], cls=cls, dsun=dsun, dist=dist, d=v,
                     lam=math.degrees(math.atan2(v[1], v[0])) % 360, bet=math.degrees(math.asin(v[2]/dist)), ra=s['ra'], de=s['de']))
rows.sort(key=lambda r: r['dist'])
pickle.dump(rows, open(os.path.join(here, 'brown_dwarfs_frame.pkl'), 'wb'))
with open(os.path.join(here, 'brown_dwarfs_in_sirius_orbit_frame.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f); w.writerow(['Name', 'Spectral type', 'RA (deg)', 'Dec (deg)', 'Distance from Sun (ly)', 'Distance from Sirius (ly)', 'X (ly)', 'Y (ly)', 'Z (ly)', 'Orbit-frame longitude (deg)', 'Orbit-frame latitude (deg)'])
    for r in rows: w.writerow([r['name'], r['spt'], r['ra'], r['de'], round(r['dsun'], 1), round(r['dist'], 1), *[round(x, 1) for x in r['d']], round(r['lam'], 1), round(r['bet'], 1)])
from collections import Counter
print(len(raw), 'entries ->', len(rows), 'systems; classes', Counter(r['cls'] for r in rows), '; within 20 ly of Sirius:', sum(r['dist'] <= 20 for r in rows))
for r in rows[:16]: print('%-18s %-10s Sun %5.1f Sirius %5.1f lon %5.1f lat %5.1f' % (r['name'], r['spt'], r['dsun'], r['dist'], r['lam'], r['bet']))
