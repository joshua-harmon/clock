# The Ursa Major moving group ("Sirius supercluster") in coordinates tied to the orbit of Sirius A and B:
# origin at the Sirius AB centre of mass, Z along the orbit's pole (the direction from which Sirius B is seen moving
# counter-clockwise), X toward Sirius B's periastron, Y = Z x X (the direction B moves at periastron).
# Members: Wikipedia list (core after King et al. 2003, plus stream stars); positions Hipparcos 2007.
import csv, math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from binary_stars import Binary, T

here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
S = Binary('Sirius')
# ---- the orbit frame, in galactic xyz
t0 = S.T0; r = S.rel3(t0); v = (S.rel3(t0 + 1e-3) - S.rel3(t0 - 1e-3))/2e-3
Zh = np.cross(r, v); Zh /= np.linalg.norm(Zh); Xh = r/np.linalg.norm(r); Yh = np.cross(Zh, Xh)
M = np.vstack([Xh, Yh, Zh])                       # galactic -> orbit frame
def lb(vec): return math.degrees(math.atan2(vec[1], vec[0])) % 360, math.degrees(math.asin(vec[2]/np.linalg.norm(vec)))
Teq = T.T                                           # galactic -> equatorial
def radec(vec): e = Teq @ vec; return math.degrees(math.atan2(e[1], e[0])) % 360, math.degrees(math.asin(e[2]/np.linalg.norm(e)))
pole_l, pole_b = lb(Zh); peri_l, peri_b = lb(Xh)
print('orbit pole: l %.1f b %.1f | RA %.1f Dec %.1f' % (pole_l, pole_b, *radec(Zh)))
print('periastron direction: l %.1f b %.1f' % (peri_l, peri_b))
print('angle between orbit plane and galactic plane: %.1f deg' % math.degrees(math.acos(abs(Zh[2]))))
LOS = S.r0                                          # Sun -> Sirius, equatorial
los_g = T @ LOS
print('tilt to our line of sight (pole vs line of sight): %.1f' % math.degrees(math.acos(abs(np.dot(Zh, los_g)))))

NAMES = {112185: ('Alioth', 'core'), 116656: ('Mizar', 'core'), 95418: ('Merak', 'core'), 103287: ('Phecda', 'core'), 106591: ('Megrez', 'core'),
         116842: ('Alcor', 'core'), 113139: ('78 UMa', 'core'), 91480: ('37 UMa', 'core'), 115043: ('HD 115043', 'core'), 109011: ('NO UMa', 'core'),
         110463: ('NP UMa', 'core'), 109647: ('DO CVn', 'core'),
         40183: ('Menkalinan (β Aur)', 'stream'), 139006: ('Alphecca (α CrB)', 'stream'), 216627: ('Skat (δ Aqr)', 'stream'), 89025: ('Adhafera (ζ Leo)', 'stream'),
         38393: ('γ Lep', 'stream'), 141003: ('β Ser', 'stream'), 129246: ('ζ Boo', 'stream'), 39587: ('χ¹ Ori', 'stream'), 87696: ('21 LMi', 'stream'),
         11171: ('χ Cet A', 'stream'), 199951: ('γ Mic', 'stream'), 102070: ('ζ Crt', 'stream'), 147584: ('ζ TrA', 'stream'), 177196: ('16 Lyr', 'stream'),
         27820: ('66 Tau', 'stream'), 180777: ('59 Dra', 'stream'), 7804: ('89 Psc', 'stream'), 75605: ('HD 75605', 'stream'), 27022: ('HD 27022', 'stream'),
         125451: ('18 Boo', 'stream'), 109799: ('HD 109799', 'stream'), 72905: ('π¹ UMa', 'stream'), 220096: ('HD 220096', 'stream'), 111397: ('29 Com', 'stream'),
         18778: ('HD 18778', 'stream'), 165185: ('HD 165185', 'stream'), 85364: ('6 Sex', 'stream'), 171746: ('HD 171746', 'stream'), 129798: ('DL Dra', 'stream'),
         112097: ('41 Vir', 'stream'), 11131: ('χ Cet B (EZ Cet)', 'stream'), 71974: ('HD 71974', 'stream'), 59747: ('DX Lyn', 'stream'), 28495: ('MS Cam', 'stream'),
         173950: ('V595 Lyr', 'stream'), None: ('Gliese 519', 'stream'), 95650: ('DS Leo', 'stream'), 238224: ('HD 238224', 'stream'), 13959: ('HD 13959', 'stream'),
         156498: ('V2369 Oph', 'stream'), 48915: ('Sirius', 'sirius')}
def xyz(ra, de, plx):
    a, d = math.radians(ra), math.radians(de); return T @ (3261.56/plx*np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]))
rows = []
for rr in csv.DictReader(open(os.path.join(here, 'data', 'uma_group_hipparcos.csv'))):
    hd = int(rr['HD']) if rr['HD'] else None
    nm, grp = NAMES[hd]
    rows.append(dict(hd=hd, hip=int(rr['HIP']), name=nm, grp=grp, V=float(rr['Vmag']), plx=float(rr['Plx']), p=xyz(float(rr['RAdeg']), float(rr['DEdeg']), float(rr['Plx']))))
sir = next(s for s in rows if s['grp'] == 'sirius'); C = sir['p']
out = []
for s in rows:
    d = M @ (s['p'] - C); dist = np.linalg.norm(d)
    if s['grp'] == 'sirius': lam = bet = float('nan')
    else: lam = math.degrees(math.atan2(d[1], d[0])) % 360; bet = math.degrees(math.asin(d[2]/dist))
    s.update(d=d, dist=dist, lam=lam, bet=bet); out.append(s)
sun = M @ (-C); sun_lam = math.degrees(math.atan2(sun[1], sun[0])) % 360; sun_bet = math.degrees(math.asin(sun[2]/np.linalg.norm(sun)))
print('Sun from Sirius in orbit frame: lon %.1f lat %.1f' % (sun_lam, sun_bet))
core = [s for s in out if s['grp'] == 'core']; cc = np.mean([s['d'] for s in core], axis=0)
print('core centre from Sirius: X %.1f Y %.1f Z %.1f ly; lon %.1f lat %.1f dist %.1f' % (*cc, math.degrees(math.atan2(cc[1], cc[0])) % 360, math.degrees(math.asin(cc[2]/np.linalg.norm(cc))), np.linalg.norm(cc)))
with open(os.path.join(here, 'uma_group_in_sirius_orbit_frame.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['Name', 'Group', 'HD', 'HIP', 'V mag', 'Distance from Sun (ly)', 'Distance from Sirius (ly)', 'X (ly, toward periastron)', 'Y (ly, along B’s motion)', 'Z (ly, toward orbit pole)',
                'Orbit-frame longitude (deg)', 'Orbit-frame latitude (deg)'])
    for s in sorted(out, key=lambda s: (s['grp'] != 'sirius', s['grp'] != 'core', s['dist'])):
        w.writerow([s['name'], {'core': 'core (nucleus)', 'stream': 'stream', 'sirius': 'reference (Sirius AB)'}[s['grp']], s['hd'] or '', s['hip'], s['V'], round(3261.56/s['plx'], 1), round(s['dist'], 1),
                    *[round(x, 1) for x in s['d']], '' if s['grp'] == 'sirius' else round(s['lam'], 1), '' if s['grp'] == 'sirius' else round(s['bet'], 1)])
np.save(os.path.join(here, 'data', 'orbit_frame_M.npy'), M)
import pickle; pickle.dump(dict(out=out, sun=(sun_lam, sun_bet, sun), cc=cc, pole=(pole_l, pole_b), M=M), open(os.path.join(here, 'data', 'orbit_frame.pkl'), 'wb'))
print(len(out), 'stars')
