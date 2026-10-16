# The zoom series turned 90°, without motion lines: the stars within 250 and 50 light years seen edge-on
# (the galactic plane as a horizontal line), and the Ursa Major nucleus close up from above and edge-on.
# Same colours, brightness and landmark labels as the top views (make_zoom_maps.py).
import ast, csv, math, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
COL = {'O': '#9bb0ff', 'B': '#aabfff', 'A': '#e6ecff', 'F': '#fff4d6', 'G': '#ffe680', 'K': '#ffb056', 'M': '#ff6a4d', 'D': '#ffffff'}
FIXCLS = {'Capella': 'G8III'}
# landmark names shared with the top views
src = open(os.path.join(here, 'make_zoom_maps.py'), encoding='utf-8').read()
LM250 = next(ast.literal_eval(n.value) for n in ast.parse(src).body if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', '') == 'LM250')

def load(fname, name_col, prop_col, cls_col, dist_col, vcol, extra=()):
    out = []
    for r in csv.DictReader(open(os.path.join(here, fname), encoding='utf-8')):
        try: l, b, d, V = float(r['Galactic l']), float(r['Galactic b']), float(r[dist_col]), float(r[vcol])
        except (ValueError, KeyError): continue
        lr, br = math.radians(l), math.radians(b)
        s = dict(name=r[name_col], proper=r.get(prop_col, ''), cls=r[cls_col], d=d, V=V,
                 X=d*math.cos(br)*math.cos(lr), Y=d*math.cos(br)*math.sin(lr), Z=d*math.sin(br))
        for k in extra: s[k] = r.get(k, '')
        out.append(s)
    return out

def dot(ax, s, scale=1.0):
    c = COL.get(FIXCLS.get(s['proper'], s['cls'])[:1].upper(), '#dddddd')
    ax.scatter([s['px']], [s['py']], s=max(4, min(160, (7 - max(s['V'], -1.5))**2*2.4))*scale, color=c, edgecolors='none', zorder=3)

def label(ax, s, text, off):
    ox, oy = off.get(text, (6, 3))
    ax.annotate(text, (s['px'], s['py']), xytext=(ox, oy), textcoords='offset points', color='#f2f2f2', fontsize=9.5, zorder=8,
                ha='right' if ox < 0 else 'left', bbox=dict(boxstyle='round,pad=.1', fc='#05070d', ec='none', alpha=.55))

def dark(fig, rect):
    ax = fig.add_axes(rect); ax.set_facecolor('#05070d')
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.tick_params(colors='#9aa', labelsize=8.5); return ax

def side_view(rows, M, rings, lab, off, fname, title, note):
    fig = plt.figure(figsize=(16, 10.5), dpi=130); fig.patch.set_facecolor('#05070d')
    ax = dark(fig, [0.06, 0.1, 0.9, 0.82])
    for s in rows: s['px'], s['py'] = s['Y'], s['Z']        # right = toward l = 90°, up = north galactic pole
    ax.axhline(0, color='#d9c27a', lw=1.2, alpha=.8, zorder=1)
    ax.text(-M*.98, M*.012, 'galactic plane, seen edge-on', color='#d9c27a', fontsize=9.5, zorder=7)
    for rr in rings:
        ax.add_patch(Circle((0, 0), rr, fill=False, ec='#d9c27a', lw=.6, ls=(0, (3, 4)), alpha=.45, zorder=1))
        ax.text(rr*.71 + M*.005, rr*.71, '%g ly' % rr, color='#cdb874', fontsize=8.5, zorder=7, clip_on=True)
    for s in sorted(rows, key=lambda s: -s['V']): dot(ax, s)
    for s in rows:
        t = lab(s)
        if t: label(ax, s, t, off)
    ax.plot(0, 0, '*', ms=16, color='#ffe24a', zorder=9); ax.annotate('Sun', (0, 0), xytext=(8, -18), textcoords='offset points', color='#ffe24a', fontsize=12, weight='bold', zorder=9)
    ax.text(-M*.98, M*.67, '↑  north galactic pole (Coma Berenices)', color='#e8e8e8', fontsize=11, ha='left', zorder=7)
    ax.text(-M*.98, -M*.69, '↓  south galactic pole (Sculptor)', color='#e8e8e8', fontsize=11, ha='left', va='bottom', zorder=7)
    ax.text(M*.98, -M*.05, 'l = 90°, Sun’s orbit (Cygnus) →', color='#e8e8e8', fontsize=10, ha='right', zorder=7)
    ax.text(-M*.98, -M*.05, '← l = 270° (Vela)', color='#e8e8e8', fontsize=10, ha='left', zorder=7)
    ax.set_xlim(-M, M); ax.set_ylim(-M*.7, M*.7); ax.set_aspect('equal')
    ax.set_xlabel('light years toward l = 90°', color='#ccc'); ax.set_ylabel('light years north of the galactic plane', color='#ccc')
    ax.set_title(title, color='white', fontsize=14.5)
    fig.text(.06, .02, note, color='#9aa', fontsize=9, wrap=True)
    fig.savefig(os.path.join(here, fname), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)

NOTE = ('The top-view maps turned 90° about the line toward the galactic centre: we look from the galactic anticentre (Auriga) toward the centre (Sagittarius), '
        'so stars are placed by how far they lie toward Cygnus or Vela (across) and how far above or below the plane of the Milky Way (up and down). '
        'Positions: Hipparcos 2007. Colours by spectral type, dot size by brightness.')

big = [s for s in load('naked_eye_stars_within_250ly.csv', 'Name', 'Proper name', 'Spectral type', 'Distance (ly)', 'V mag', ('Ursa Major group nucleus',)) if s['d'] <= 250]
def lab250(s):
    if s['d'] < 50: return ''
    if s['proper'] == 'Merak': return 'Big Dipper core'
    return LM250.get(s['name'], '')
OFF250 = {'Big Dipper core': (6, 3), 'Dubhe': (-6, 4), 'Alkaid': (6, -10), 'Aldebaran': (-6, -10), 'Castor': (6, 3), 'Regulus': (-6, 4)}
side_view(big, 260, [50, 100, 150, 200, 250], lab250, OFF250, 'side_view_within_250ly.png',
          'Naked-eye stars within 250 light years (%d), turned 90°: seen edge-on' % len(big),
          NOTE + ' Stars nearer than 50 ly are labelled in the 50 ly view.')

fifty = load('stars_within_50ly.csv', 'Catalog Name', 'Proper Name', 'Classification', 'Dist (ly)', 'Visual Mag')
def lab50(s):
    p = s['proper']
    if s['name'] == 'Zeta Reticuli': return 'Zeta Reticuli'
    if p in ('Sol', ''): return s['name'] if s['V'] < 3.0 else ''
    return 'Alpha Centauri' if p == 'Rigil Kentaurus' else p
OFF50 = {'Sirius': (-6, -12), 'Procyon': (-6, 4), 'Alpha Centauri': (-4, 6), 'Altair': (6, -12), 'Arcturus': (6, 3), 'Muphrid': (-6, 4),
         'Denebola': (-6, 4), 'Capella': (6, 3), 'Fomalhaut': (6, -12), 'Vega': (6, 3), 'Keid': (-6, -10), 'Talita': (-6, 4), 'Achird': (6, 3), 'Errai': (6, 3)}
side_view(fifty, 52, [10, 20, 30, 40, 50], lab50, OFF50, 'side_view_within_50ly.png',
          'Naked-eye stars within 50 light years (%d systems), turned 90°: seen edge-on' % len(fifty), NOTE)

# ---- the Ursa Major nucleus close up, with its neighbours, from above and edge-on
core = [s for s in big if s.get('Ursa Major group nucleus') == 'yes']
cx, cy, cz = (sum(s[k] for s in core)/len(core) for k in ('X', 'Y', 'Z'))
R = 28
box = [s for s in big if abs(s['X'] - cx) < R and abs(s['Y'] - cy) < R and abs(s['Z'] - cz) < R]
fig = plt.figure(figsize=(18, 9.6), dpi=130); fig.patch.set_facecolor('#05070d')
for k, (rect, horiz, vert, xl, yl, ttl) in enumerate([
        ([0.04, 0.1, 0.44, 0.78], 'Y', 'X', 'light years toward l = 90° (Cygnus)', 'light years toward the galactic centre (down)', 'Seen from above the galactic plane'),
        ([0.54, 0.1, 0.44, 0.78], 'Y', 'Z', 'light years toward l = 90° (Cygnus)', 'light years north of the galactic plane', 'Turned 90°: seen edge-on')]):
    ax = dark(fig, rect)
    c0 = {'X': cx, 'Y': cy, 'Z': cz}
    for s in box:
        s['px'] = s[horiz] - c0[horiz]
        s['py'] = -(s[vert] - c0[vert]) if vert == 'X' else s[vert] - c0[vert]
    for s in sorted(box, key=lambda s: -s['V']): dot(ax, s, 1.7)
    for s in core:
        ax.add_patch(Circle((s['px'], s['py']), 1.2, fill=False, ec='#d0302a', lw=1.2, zorder=4))
    for s in box:
        t = s['proper'] or (s['name'] if s['V'] < 4.0 or s.get('Ursa Major group nucleus') == 'yes' else '')
        if t: label(ax, s, t, {'Alcor': (-6, -12), 'Mizar': (6, 3), '78 UMa': (6, -12)})
    ax.plot(0, 0, '+', ms=14, color='#9fb6d8', zorder=5)
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect('equal'); ax.grid(color='#1c2636', lw=.6)
    ax.set_xlabel(xl, color='#ccc'); ax.set_ylabel(yl, color='#ccc'); ax.set_title(ttl, color='white', fontsize=12.5)
fig.suptitle('The Ursa Major group nucleus close up: %d naked-eye stars in a %d-light-year box about 80 ly from the Sun (red rings: the nucleus)' % (len(box), 2*R),
             color='white', fontsize=13.5)
fig.text(.04, .02, 'Centred on the middle of the nucleus (+). The seven nucleus stars (Merak, Phecda, Megrez, Alioth, Mizar, Alcor, 78 UMa) share one motion through space; '
         'Alkaid, at the handle end of the Big Dipper, falls inside the box but moves differently and is not a member; Dubhe, at the other end, lies outside the box at 123 ly. Hipparcos 2007.', color='#9aa', fontsize=9, wrap=True)
fig.savefig(os.path.join(here, 'ursa_major_nucleus_close_up.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)
print('done', len(big), len(fifty), len(box), [s['proper'] or s['name'] for s in box])
