# Top views (from the north galactic pole) of the stars within 250 and 50 light years, in the same style as the
# Milky Way and Orion Arm maps: galactic longitude drawn from the Sun, distance rings, landmark stars labelled.
# Positions are the real Hipparcos positions projected onto the galactic plane.
import csv, math, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
COL = {'O': '#9bb0ff', 'B': '#aabfff', 'A': '#e6ecff', 'F': '#fff4d6', 'G': '#ffe680', 'K': '#ffb056', 'M': '#ff6a4d', 'D': '#ffffff'}

def load(fname, name_col, prop_col, cls_col, dist_col, vcol):
    out = []
    for r in csv.DictReader(open(os.path.join(here, fname), encoding='utf-8')):
        try: l, b, d, V = float(r['Galactic l']), float(r['Galactic b']), float(r[dist_col]), float(r[vcol])
        except (ValueError, KeyError): continue
        lr, br = math.radians(l), math.radians(b); dp = d*math.cos(br)
        out.append(dict(name=r[name_col], proper=r.get(prop_col, ''), cls=r[cls_col], d=d, V=V, z=d*math.sin(br),
                        x=dp*math.sin(lr), y=-dp*math.cos(lr)))     # x toward l = 90°, y away from the galactic centre
    return out

def zoom_map(rows, M, rings, label, offsets, fname, title, note, extra=None):
    fig = plt.figure(figsize=(15, 14), dpi=130); fig.patch.set_facecolor('#05070d')
    ax = fig.add_axes([0.03, 0.07, 0.94, 0.88]); ax.set_facecolor('#05070d')
    for l in range(0, 360, 15):
        lr = math.radians(l); L = M*1.6; major = l % 30 == 0
        ax.plot([0, L*math.sin(lr)], [0, -L*math.cos(lr)], color='#d9c27a', lw=1.1 if l % 90 == 0 else (.7 if major else .35), alpha=.75 if major else .4, zorder=1)
        if major:
            s, c = math.sin(lr), -math.cos(lr); t = min(M*.965/abs(s) if abs(s) > 1e-9 else 1e9, M*.965/abs(c) if abs(c) > 1e-9 else 1e9)
            ax.text(t*s, t*c, '%d°' % l, color='#f0d98a', fontsize=11 if l % 90 == 0 else 9.5, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=.2', fc='#05070d', ec='none', alpha=.85), zorder=7)
    for rr in rings:
        ax.add_patch(Circle((0, 0), rr, fill=False, ec='#d9c27a', lw=.6, ls=(0, (3, 4)), alpha=.55, zorder=1))
        ax.text(-rr*.97, rr*.26, '%g ly' % rr, color='#cdb874', fontsize=8.5, zorder=7)
    for r in sorted(rows, key=lambda r: -r['V']):
        c = COL.get(r['cls'][:1].upper(), '#dddddd'); s = max(4, min(150, (7 - max(r['V'], -1.5))**2*2.4))
        ax.scatter([r['x']], [r['y']], s=s, color=c, edgecolors='none', zorder=3, alpha=.95)
        lab = label(r)
        if lab:
            ox, oy = offsets.get(lab, (6, 3))
            ax.annotate(lab, (r['x'], r['y']), xytext=(ox, oy), textcoords='offset points', color='#f2f2f2', fontsize=9.5, zorder=8,
                        ha='right' if ox < 0 else 'left', bbox=dict(boxstyle='round,pad=.1', fc='#05070d', ec='none', alpha=.55))
    if extra: extra(ax)
    ax.plot(0, 0, marker='*', ms=16, color='#ffe24a', zorder=9)
    ax.annotate('Sun', (0, 0), xytext=(8, -15), textcoords='offset points', color='#ffe24a', fontsize=12, weight='bold', zorder=9)
    ax.text(0, -M*.9, '↓  toward the galactic centre (Sagittarius), 26,600 ly', color='#e8e8e8', fontsize=11, ha='center', zorder=7)
    ax.text(0, M*.88, 'away from the galactic centre (Auriga)  ↑', color='#e8e8e8', fontsize=11, ha='center', zorder=7)
    ax.text(M*.9, M*.05, 'Sun’s orbit →\n(Cygnus)', color='#e8e8e8', fontsize=10, ha='right', zorder=7)
    ax.text(-M*.9, M*.05, '(Vela)', color='#e8e8e8', fontsize=10, ha='left', zorder=7)
    for k, (t, key) in enumerate([('O/B', 'B'), ('A', 'A'), ('F', 'F'), ('G', 'G'), ('K', 'K'), ('M', 'M')]):
        ax.scatter([-M*.95 + k*M*.11], [-M*.8], s=45, color=COL[key]); ax.text(-M*.93 + k*M*.11, -M*.8, t, color='#ccc', fontsize=9, va='center')
    ax.text(-M*.95, -M*.85, 'star colours by spectral type; dot size by brightness', color='#aaa', fontsize=8.5)
    ax.set_xlim(-M, M); ax.set_ylim(-M, M); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, color='white', fontsize=15)
    fig.text(.03, .03, note, color='#9aa', fontsize=9, wrap=True)
    fig.savefig(os.path.join(here, fname), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)

NOTE = ('Seen from the north galactic pole; positions projected onto the galactic plane (Hipparcos 2007). Galactic longitude: 0° toward the galactic centre, '
        '90° in the direction of the Sun’s orbit, 180° away from the centre, 270° toward Vela. All of this lies deep inside the Local (Orion) Arm, which is about 3,000 ly wide.')

# ---- 250 light years: every naked-eye star, landmarks labelled
big = [r for r in load('naked_eye_stars_within_250ly.csv', 'Name', 'Proper name', 'Spectral type', 'Distance (ly)', 'V mag') if r['d'] <= 250]
LM250 = {'Alpha CMa': 'Sirius', 'Alpha Boo': 'Arcturus', 'Alpha Lyr': 'Vega', 'Alpha Aur': 'Capella', 'Alpha Tau': 'Aldebaran', 'Alpha Eri': 'Achernar',
         'Alpha Vir': 'Spica', 'Alpha Leo': 'Regulus', 'Alpha PsA': 'Fomalhaut', 'Alpha Aql': 'Altair', 'Alpha CMi': 'Procyon', 'Beta Gem': 'Pollux',
         'Alpha Gem': 'Castor', 'Alpha UMa': 'Dubhe', 'Eta UMa': 'Alkaid', 'Alpha Car': 'Canopus', 'Alpha Cen': 'Alpha Centauri', 'Alpha Oph': 'Rasalhague',
         'Alpha Ari': 'Hamal', 'Alpha CrB': 'Alphecca', 'Alpha Hya': 'Alphard', 'Gamma Dra': 'Eltanin', 'Alpha Peg': 'Markab', 'Beta Leo': 'Denebola',
         'Alpha Gru': 'Alnair', 'Alpha Pav': 'Peacock', 'Beta UMi': 'Kochab', 'Alpha Cep': 'Alderamin', 'Alpha And': 'Alpheratz', 'Gamma Leo': 'Algieba',
         'Beta Per': 'Algol', 'Alpha Tri': 'Mothallah', 'Alpha Ser': 'Unukalhai'}
def lab250(r):
    if r['proper'] == 'Merak': return 'Big Dipper core'
    if r['d'] < 50: return ''          # the nearest stars are labelled in the 50 ly map
    return LM250.get(r['name'], '') if not (r['name'] == 'Alpha Cen' and r['V'] > 1) else ''
OFF250 = {'Sirius': (6, -10), 'Procyon': (6, 3), 'Alpha Centauri': (-6, 4), 'Altair': (6, -10), 'Vega': (6, 4), 'Fomalhaut': (-6, -10),
          'Pollux': (-6, 4), 'Castor': (6, 3), 'Arcturus': (6, 3), 'Capella': (-6, 4), 'Aldebaran': (6, -10), 'Big Dipper core': (6, 3),
          'Dubhe': (-6, 4), 'Alkaid': (6, 3), 'Denebola': (6, -10), 'Rasalhague': (6, 3)}
def hyades(ax):
    ax.add_patch(Circle((0, 0), 50, fill=False, ec='#9fb6d8', lw=1.2, zorder=6))
    ax.text(0, 52, '50 ly: see the next map (Sirius, Vega, Arcturus, Altair, Procyon, Fomalhaut, Capella, Pollux …)', color='#9fb6d8', fontsize=8.5, ha='center', va='bottom', zorder=8,
            bbox=dict(boxstyle='round,pad=.15', fc='#05070d', ec='none', alpha=.7))
    lr = math.radians(180.1); x, y = 153*math.sin(lr), -153*math.cos(lr)
    ax.add_patch(Circle((x, y), 14, fill=False, ec='#9fe0ff', lw=1, zorder=6)); ax.annotate('Hyades cluster', (x, y), xytext=(-12, 10), textcoords='offset points',
                                                                                    color='#9fe0ff', fontsize=9.5, ha='right', zorder=8)
zoom_map(big, 260, [50, 100, 150, 200, 250], lab250, OFF250, 'top_view_within_250ly.png',
         'Naked-eye stars within 250 light years (%d), seen from above the galactic plane' % len(big), NOTE, hyades)

# ---- 50 light years: every naked-eye system, named stars labelled
fifty = load('stars_within_50ly.csv', 'Catalog Name', 'Proper Name', 'Classification', 'Dist (ly)', 'Visual Mag')
def lab50(r):
    p = r['proper']
    if r['name'] == 'Zeta Reticuli': return 'Zeta Reticuli'
    if p in ('Sol', ''): return r['name'] if r['V'] < 3.0 else ''
    return 'Alpha Centauri' if p == 'Rigil Kentaurus' else p
OFF50 = {'Sirius': (6, -10), 'Procyon': (6, 3), 'Alpha Centauri': (-6, 4), 'Altair': (6, 3), 'Alshain': (6, -10), 'Achird': (-6, 4),
         'Fomalhaut': (6, 3), 'Keid': (-6, 4), 'Pollux': (-6, 4), 'Capella': (6, 3), 'Talita': (6, 3), 'Al Haud': (-6, 4), 'Zavijava': (6, -10),
         'Denebola': (-6, 4), 'Porrima': (6, 3), 'Muphrid': (6, 3), 'Arcturus': (-6, -12), 'Alsafi': (6, 3), 'Errai': (-6, 4)}
zoom_map(fifty, 52, [10, 20, 30, 40, 50], lab50, OFF50, 'top_view_within_50ly.png',
         'Naked-eye stars within 50 light years (%d systems), seen from above the galactic plane' % len(fifty), NOTE)
print('done', len(big), len(fifty))
