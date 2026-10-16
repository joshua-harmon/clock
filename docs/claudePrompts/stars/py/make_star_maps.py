# Builds the nearby-star list and maps from hipparcos_naked_eye_within_250ly.csv
# (Hipparcos 2007 reduction I/311 + Hipparcos main catalogue I/239 + HD/HR/Bayer/Flamsteed index IV/27A, via VizieR).
# Usage: python make_star_maps.py path\to\hipparcos_naked_eye_within_250ly.csv
import csv, math, sys, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else 'hipparcos_naked_eye_within_250ly.csv'
out = os.path.dirname(os.path.abspath(__file__))
GREEK = dict(alf='Alpha', bet='Beta', gam='Gamma', del_='Delta', eps='Epsilon', zet='Zeta', eta='Eta', the='Theta', iot='Iota',
             kap='Kappa', lam='Lambda', mu='Mu', nu='Nu', ksi='Xi', omi='Omicron', pi='Pi', rho='Rho', sig='Sigma', tau='Tau',
             ups='Upsilon', phi='Phi', chi='Chi', psi='Psi', ome='Omega')
PROPER = {32349:'Sirius', 91262:'Vega', 97649:'Altair', 69673:'Arcturus', 37279:'Procyon', 24608:'Capella', 37826:'Pollux',
          113368:'Fomalhaut', 71683:'Rigil Kentaurus', 71681:'Toliman', 53910:'Merak', 58001:'Phecda', 59774:'Megrez',
          62956:'Alioth', 65378:'Mizar', 65477:'Alcor', 54061:'Dubhe', 67301:'Alkaid', 49669:'Regulus', 57632:'Denebola',
          36850:'Castor', 21421:'Aldebaran', 7588:'Achernar', 65474:'Spica', 107556:'Deneb Algedi', 105199:'Alderamin', 86032:'Rasalhague'}
# Ursa Major group nucleus (King et al. 2003, AJ 125, 1980): the five middle Big Dipper stars, Alcor and 78 UMa
UMA_CORE = {('bet','UMa'), ('gam','UMa'), ('del','UMa'), ('eps','UMa'), ('zet','UMa'), ('80','UMa'), ('78','UMa')}

def gal(ra, dec):
    # ICRS (J2000) to galactic, standard IAU rotation
    T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
                  [ 0.4941094279, -0.4448296300,  0.7469822445],
                  [-0.8676661490, -0.1980763734,  0.4559837762]])
    a, d = math.radians(ra), math.radians(dec)
    v = T @ np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
    l = math.degrees(math.atan2(v[1], v[0])) % 360
    b = math.degrees(math.asin(max(-1, min(1, v[2]))))
    return l, b

rows, seen = [], set()
for r in csv.DictReader(open(src, encoding='utf-8')):
    h = int(r['HIP'])
    if h in seen: continue
    seen.add(h)
    ra, dec, plx = float(r['RArad']), float(r['DErad']), float(r['Plx'])
    l, b = gal(ra, dec)
    bay, fl, cst = r['Bayer'].strip(), r['Fl'].strip(), r['Cst'].strip()
    g = GREEK.get('del_' if bay[:3] == 'del' else bay.rstrip('0123456789'), None)
    if g and cst: bname = g + (bay[-1] if bay[-1].isdigit() else '') + ' ' + cst
    elif bay and cst: bname = bay + ' ' + cst
    elif fl and cst: bname = fl + ' ' + cst
    elif r['HR'].strip(): bname = 'HR ' + r['HR'].strip()
    else: bname = 'HIP %d' % h
    core = (bay[:3], cst) in UMA_CORE or (fl, cst) in UMA_CORE
    rows.append(dict(HIP=h, name=bname, proper=PROPER.get(h, ''), ra=ra, dec=dec, l=l, b=b, V=float(r['Vmag']),
                     sp=r['SpType'].strip(), plx=plx, eplx=float(r['e_Plx']), ly=3261.56/plx,
                     pmra=float(r['pmRA']), pmde=float(r['pmDE']), uma=core))
rows.sort(key=lambda r: r['ly'])

def hms(ra):
    h = ra/15; return '%02d %04.1f' % (int(h), (h-int(h))*60)
def dms(d):
    s = '-' if d < 0 else '+'; d = abs(d); return '%s%02d %02d' % (s, int(d), int(round((d-int(d))*60)) % 60)

with open(os.path.join(out, 'naked_eye_stars_within_250ly.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['HIP','Name','Proper name','RA (h m)','Dec (d m)','Galactic l','Galactic b','V mag','Spectral type',
                'Parallax (mas)','Parallax error','Distance (ly)','pm RA (mas/yr)','pm Dec (mas/yr)','Ursa Major group nucleus'])
    for r in rows:
        w.writerow([r['HIP'], r['name'], r['proper'], hms(r['ra']), dms(r['dec']), '%.1f' % r['l'], '%+.1f' % r['b'], r['V'], r['sp'],
                    r['plx'], r['eplx'], '%.1f' % r['ly'], r['pmra'], r['pmde'], 'yes' if r['uma'] else ''])

n100 = sum(r['ly'] <= 100 for r in rows)
print('stars', len(rows), 'within 100 ly', n100, 'UMa core', [r['name'] for r in rows if r['uma']])

# ---- figure 1: galactic sky map (Aitoff), centre of the Galaxy in the middle, longitude increasing to the left
def sky(maxly, fname, title):
    S = [r for r in rows if r['ly'] <= maxly]
    fig = plt.figure(figsize=(14, 7.6), dpi=130); ax = fig.add_subplot(111, projection='aitoff')
    L = np.array([-(((r['l'] + 180) % 360) - 180) for r in S]); B = np.array([r['b'] for r in S])
    sz = np.array([max(2, (6.6 - r['V'])**2 * 3.2) for r in S]); dist = np.array([r['ly'] for r in S])
    sc = ax.scatter(np.radians(L), np.radians(B), s=sz, c=dist, cmap='viridis_r', vmin=0, vmax=maxly, edgecolors='none', alpha=.9)
    U = [r for r in S if r['uma']]
    ax.scatter([math.radians(-(((r['l']+180)%360)-180)) for r in U], [math.radians(r['b']) for r in U], s=150, facecolors='none', edgecolors='#d0302a', linewidths=1.6, label='Ursa Major group nucleus')
    for r in S:
        if r['proper'] and r['proper'] not in ('Toliman','Phecda','Megrez','Alioth','Alcor','Castor') and (r['V'] < 1.6 or r['proper'] in ('Merak','Mizar','Sirius','Dubhe','Alkaid')):
            x = math.radians(-(((r['l']+180)%360)-180)); y = math.radians(r['b'])
            ax.annotate({'Merak':'Merak  (Big Dipper core)','Pollux':'Pollux, Castor','Rigil Kentaurus':'Alpha Centauri'}.get(r['proper'],r['proper']), (x, y), xytext=(5, 4) if r['proper']!='Merak' else (-8,10), textcoords='offset points', fontsize=8, color='#222', ha='right' if r['proper']=='Merak' else 'left')
    ax.set_xticklabels(['150°','120°','90°','60°','30°','0°','330°','300°','270°','240°','210°'], fontsize=8)
    ax.grid(True, color='#bbb', lw=.5); cb = fig.colorbar(sc, ax=ax, orientation='horizontal', fraction=.04, pad=.07, shrink=.5)
    cb.set_label('distance (light years)')
    ax.set_title(title, fontsize=13, pad=14); ax.legend(loc='lower right', fontsize=8, frameon=False)
    fig.text(.5, .015, 'Galactic coordinates: galactic centre in the middle, longitude increasing to the left, plane of the Milky Way along the equator. Dot size: brightness. Hipparcos (van Leeuwen 2007).', ha='center', fontsize=8, color='#555')
    fig.savefig(os.path.join(out, fname), bbox_inches='tight'); plt.close(fig)

sky(100, 'sky_galactic_within_100ly.png', 'Naked-eye stars within 100 light years (%d), galactic coordinates' % n100)
sky(250, 'sky_galactic_within_250ly.png', 'Naked-eye stars within 250 light years (%d), galactic coordinates' % len(rows))

# ---- figure 2: positions around the Sun, top view (galactic plane) and side view
def xyz(r):
    lr, br = math.radians(r['l']), math.radians(r['b']); d = r['ly']
    return d*math.cos(br)*math.cos(lr), d*math.cos(br)*math.sin(lr), d*math.sin(br)
P = np.array([xyz(r) for r in rows]); V = np.array([r['V'] for r in rows]); um = np.array([r['uma'] for r in rows])
fig, axs = plt.subplots(1, 2, figsize=(15, 7.6), dpi=130)
for ax, (i, j, lab) in zip(axs, [(0, 1, ('X: toward galactic centre (ly)', 'Y: direction of the Sun’s orbit (ly)')), (0, 2, ('X: toward galactic centre (ly)', 'Z: toward north galactic pole (ly)'))]):
    sz = np.clip((6.6 - V)**2 * 2.2, 1.5, None)
    ax.scatter(P[:, i], P[:, j], s=sz, c='#3b6fb0', alpha=.55, edgecolors='none')
    ax.scatter(P[um, i], P[um, j], s=70, facecolors='none', edgecolors='#d0302a', lw=1.5, label='Ursa Major group nucleus')
    for R_ in (100, 250):
        t = np.linspace(0, 2*np.pi, 300); ax.plot(R_*np.cos(t), R_*np.sin(t), color='#999', lw=.8, ls='--')
        ax.text(R_*.72, R_*.72, '%d ly' % R_, color='#777', fontsize=8)
    ax.plot(0, 0, marker='*', color='#e8a400', ms=13, zorder=5); ax.text(6, -14, 'Sun', fontsize=9)
    for k, r in enumerate(rows):
        if r['proper'] in ('Sirius', 'Arcturus', 'Capella', 'Merak', 'Aldebaran', 'Spica', 'Achernar', 'Dubhe', 'Alkaid'):
            ax.annotate(r['proper'], (P[k, i], P[k, j]), xytext=(4, 3), textcoords='offset points', fontsize=7.5)
    ax.set_xlim(-260, 260); ax.set_ylim(-260, 260); ax.set_aspect('equal'); ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1]); ax.grid(color='#eee')
axs[0].set_title('Top view: the galactic plane seen from the north'); axs[1].set_title('Side view: up and down out of the plane')
axs[0].legend(loc='lower left', fontsize=8, frameon=False)
fig.suptitle('Naked-eye stars within 250 light years around the Sun (%d stars)' % len(rows), fontsize=13)
fig.savefig(os.path.join(out, 'positions_around_sun_250ly.png'), bbox_inches='tight'); plt.close(fig)
print('done')
