# Well-known non-member stars (the Sun, Barnard's Star, Procyon, ...) in the Sirius A-B orbit frame.
import csv, math, os, pickle
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(here, 'orbit_frame.pkl'), 'rb')); M = D['M']
C = -np.linalg.solve(M, D['sun'][2])          # Sirius position (galactic, ly) from the stored Sun vector
WANT = {('stars_within_10ly.csv', 'Barnard’s Star'): 'Barnard’s Star', ('stars_within_10ly.csv', 'Wolf 359'): 'Wolf 359',
        ('stars_within_50ly.csv', 'Alpha Centauri'): 'Alpha Centauri', ('stars_within_50ly.csv', 'Alpha Canis Minoris'): 'Procyon',
        ('stars_within_50ly.csv', '61 Cygni'): '61 Cygni', ('stars_within_50ly.csv', 'Tau Ceti'): 'Tau Ceti', ('stars_within_50ly.csv', 'Epsilon Eridani'): 'Epsilon Eridani',
        ('stars_within_50ly.csv', 'Alpha Aquilae'): 'Altair', ('stars_within_50ly.csv', 'Alpha Lyrae'): 'Vega', ('stars_within_50ly.csv', 'Alpha Piscis Austrini'): 'Fomalhaut',
        ('stars_within_50ly.csv', 'Beta Geminorum'): 'Pollux', ('stars_within_50ly.csv', 'Alpha Boötis'): 'Arcturus', ('stars_within_50ly.csv', 'Alpha Aurigae'): 'Capella'}
refs = [dict(name='Sun', V=-26.7, p=np.zeros(3))]
for (f, cat), nm in WANT.items():
    for r in csv.DictReader(open(os.path.join(here, f), encoding='utf-8')):
        if r['Catalog Name'] == cat:
            l, b, d = math.radians(float(r['Galactic l'])), math.radians(float(r['Galactic b'])), float(r['Dist (ly)'])
            refs.append(dict(name=nm, V=float(r['Visual Mag']), p=d*np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)]))); break
    else: print('not found', nm)
rows = []
for s in refs:
    v = M @ (s['p'] - C); dist = np.linalg.norm(v)
    s.update(d=v, dist=dist, lam=math.degrees(math.atan2(v[1], v[0])) % 360, bet=math.degrees(math.asin(v[2]/dist)), sundist=np.linalg.norm(s['p']))
pickle.dump(refs, open(os.path.join(here, 'reference_stars_frame.pkl'), 'wb'))
with open(os.path.join(here, 'reference_stars_in_sirius_orbit_frame.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f); w.writerow(['Name', 'Distance from Sun (ly)', 'Distance from Sirius (ly)', 'X (ly, toward periastron)', 'Y (ly, along B’s motion)', 'Z (ly, toward orbit pole)', 'Orbit-frame longitude (deg)', 'Orbit-frame latitude (deg)'])
    for s in sorted(refs, key=lambda s: s['dist']):
        w.writerow([s['name'], round(s['sundist'], 1), round(s['dist'], 1), *[round(x, 1) for x in s['d']], round(s['lam'], 1), round(s['bet'], 1)])
        print('%-16s Sun %5.1f  Sirius %5.1f  X %6.1f Y %6.1f Z %6.1f  lon %5.1f lat %5.1f' % (s['name'], s['sundist'], s['dist'], *s['d'], s['lam'], s['bet']))
