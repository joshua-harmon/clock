# Motion maps from the 3-D velocities (U, V, W) in ursa_major_stream_velocities.csv:
#  1. top view of the stars within 50 ly with motion lines, in two reference frames side by side:
#     relative to the Sun, and relative to the average motion of the nearby stars (the Local Standard of Rest, LSR)
#  2. the same stars turned 90°: the edge-on (side) view, up and down out of the galactic plane
#  3. the Ursa Major group nucleus close up, with each star's motion relative to the group's own average
import csv, math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
COL = {'O': '#9bb0ff', 'B': '#aabfff', 'A': '#e6ecff', 'F': '#fff4d6', 'G': '#ffe680', 'K': '#ffb056', 'M': '#ff6a4d', 'D': '#ffffff'}
SUN_LSR = np.array([11.1, 12.24, 7.25])     # the Sun's motion relative to the LSR, km/s (Schönrich, Binney & Dehnen 2010)
KMS = 1.0227e-6 * 3.26156                     # 1 km/s in light years per year

stars = []
for r in csv.DictReader(open(os.path.join(here, 'ursa_major_stream_velocities.csv'), encoding='utf-8')):
    l, b, d = math.radians(float(r['Galactic l'])), math.radians(float(r['Galactic b'])), float(r['Distance (ly)'])
    stars.append(dict(HIP=int(r['HIP']), name=r['Proper name'] or r['Name'], cls=r['Spectral type'], V=float(r['V mag']), d=d,
                      X=d*math.cos(b)*math.cos(l), Y=d*math.cos(b)*math.sin(l), Z=d*math.sin(b),
                      v=np.array([float(r['U (km/s)']), float(r['V (km/s)']), float(r['W (km/s)'])]), match=r['Stream match']))
print(len(stars), 'stars with 3-D velocities')

def size(V): return max(6, min(160, (7 - max(V, -1.5))**2*2.6))
FIX = {'Capella': 'G8III'}   # Hipparcos lists the Capella pair under the companion's type
def colour(c): return COL.get(c[:1].upper(), '#dddddd')

def frame(ax, M, xlabel, ylabel):
    ax.set_facecolor('#05070d'); ax.set_xlim(-M, M); ax.set_ylim(-M, M); ax.set_aspect('equal')
    for s in ax.spines.values(): s.set_color('#3a4a60')
    ax.tick_params(colors='#9aa', labelsize=8); ax.set_xlabel(xlabel, color='#ccc'); ax.set_ylabel(ylabel, color='#ccc')
    ax.grid(color='#1c2636', lw=.6)

def draw(ax, rows, px, py, vel, years, labels=True, lab=lambda s: s['name'] if s['V'] < 2.6 or s['name'] in ('Sirius', 'Procyon', 'Fomalhaut', 'Capella') else ''):
    for s in sorted(rows, key=lambda s: -s['V']):
        x, y = px(s), py(s); dx, dy = px_v(vel(s), px), py_v(vel(s), py)
        c = colour(FIX.get(s['name'], s['cls']))
        ax.plot([x, x + dx*years*KMS], [y, y + dy*years*KMS], color=c, lw=1.0, alpha=.85, zorder=2)
        ax.scatter([x], [y], s=size(s['V']), color=c, edgecolors='none', zorder=3)
        if labels and lab(s):
            ax.annotate(lab(s), (x, y), xytext=(5, 3), textcoords='offset points', color='#eee', fontsize=8.5, zorder=4)

# which velocity component goes with which position axis
def px_v(v, f): return v[{'X': 0, 'Y': 1, 'Z': 2}[f.axis]] * f.sign
def py_v(v, f): return v[{'X': 0, 'Y': 1, 'Z': 2}[f.axis]] * f.sign
class Ax:
    def __init__(self, axis, sign): self.axis, self.sign = axis, sign
    def __call__(self, s): return s[self.axis] * self.sign
# top view, same orientation as the other maps: right = toward l = 90° (+Y), up = away from the centre (-X)
TOPX, TOPY = Ax('Y', 1), Ax('X', -1)
# side view (the map turned 90° about the line toward the galactic centre): right = l = 90° (+Y), up = north galactic pole (+Z)
SIDEX, SIDEY = Ax('Y', 1), Ax('Z', 1)

near = [s for s in stars if s['d'] <= 50]
YRS = 50000
fig, axs = plt.subplots(1, 2, figsize=(20, 10.5), dpi=120); fig.patch.set_facecolor('#05070d')
for ax, (title, vel, sunv) in zip(axs, [('Relative to the Sun', lambda s: s['v'], None),
                                        ('Relative to the average of the nearby stars (Local Standard of Rest)', lambda s: s['v'] + SUN_LSR, SUN_LSR)]):
    frame(ax, 55, 'toward l = 90°, the direction of the Sun’s orbit (ly)', 'away from the galactic centre (ly)')
    draw(ax, near, TOPX, TOPY, vel, YRS)
    ax.plot(0, 0, '*', ms=15, color='#ffe24a', zorder=5)
    if sunv is not None:
        ax.annotate('', xy=(sunv[1]*YRS*KMS, -sunv[0]*YRS*KMS), xytext=(0, 0), arrowprops=dict(arrowstyle='-|>', color='#ffe24a', lw=2), zorder=6)
        ax.text(sunv[1]*YRS*KMS + 1, -sunv[0]*YRS*KMS - 3, 'the Sun’s own motion', color='#ffe24a', fontsize=9)
    ax.set_title(title, color='white', fontsize=12.5)
fig.suptitle('Stars within 50 light years seen from above: where each moves in the next %s years (lines), %d stars with measured speeds' % (format(YRS, ','), len(near)),
             color='white', fontsize=14)
fig.text(.5, .015, 'Line = straight-line motion over %s years, in the colour of the star; dot size = brightness. Left: as seen from the moving Sun. Right: the Sun’s own motion '
         '(about 18 km/s toward Hercules and Vega) added back, so the stars move about the local average and the Sun moves too. Galactic centre is down.' % format(YRS, ','),
         color='#9aa', fontsize=9, ha='center', wrap=True)
fig.savefig(os.path.join(here, 'motions_within_50ly_two_frames.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)

# --- side view, turned 90°
fig, axs = plt.subplots(1, 2, figsize=(20, 10.5), dpi=120); fig.patch.set_facecolor('#05070d')
for ax, (title, vel, sunv) in zip(axs, [('Relative to the Sun', lambda s: s['v'], None),
                                        ('Relative to the Local Standard of Rest', lambda s: s['v'] + SUN_LSR, SUN_LSR)]):
    frame(ax, 55, 'toward l = 90°, the direction of the Sun’s orbit (ly)', 'north of the galactic plane (ly)')
    ax.axhline(0, color='#d9c27a', lw=.8, alpha=.6); ax.text(-53, 1, 'galactic plane (seen edge-on)', color='#d9c27a', fontsize=8.5)
    draw(ax, near, SIDEX, SIDEY, vel, YRS)
    ax.plot(0, 0, '*', ms=15, color='#ffe24a', zorder=5)
    if sunv is not None:
        ax.annotate('', xy=(sunv[1]*YRS*KMS, sunv[2]*YRS*KMS), xytext=(0, 0), arrowprops=dict(arrowstyle='-|>', color='#ffe24a', lw=2), zorder=6)
        ax.text(sunv[1]*YRS*KMS + 1, sunv[2]*YRS*KMS, 'the Sun’s own motion', color='#ffe24a', fontsize=9)
    ax.set_title(title, color='white', fontsize=12.5)
fig.suptitle('The same stars turned 90°: the side view, looking from the galactic anticentre toward the galactic centre', color='white', fontsize=14)
fig.text(.5, .015, 'Up = toward the north galactic pole (Coma Berenices); right = direction of the Sun’s orbit (Cygnus). Lines: motion over %s years. '
         'Relative to the local average, the Sun is rising out of the plane at about 7 km/s; it crosses the plane every 30–35 million years.' % format(YRS, ','),
         color='#9aa', fontsize=9, ha='center', wrap=True)
fig.savefig(os.path.join(here, 'motions_within_50ly_side_view.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)

# --- the Ursa Major group nucleus in its own frame
core = [s for s in stars if s['match'] == 'nucleus']
c0 = np.mean([[s['X'], s['Y'], s['Z']] for s in core], axis=0); v0 = np.mean([s['v'] for s in core], axis=0)
for s in core: s['rx'], s['ry'], s['rz'] = s['X'] - c0[0], s['Y'] - c0[1], s['Z'] - c0[2]
spread = max(abs(v) for s in core for v in (s['rx'], s['ry'], s['rz']))
print('nucleus centre %.1f ly from the Sun; members within %.1f ly of it; mean velocity %s; speeds about the mean: %s' %
      (np.linalg.norm(c0), spread, np.round(v0, 1), [round(float(np.linalg.norm(s['v'] - v0)), 1) for s in core]))
MY = 1_000_000
fig, axs = plt.subplots(1, 2, figsize=(18, 9.5), dpi=120); fig.patch.set_facecolor('#05070d')
for ax, (a, b, xl, yl, ttl) in zip(axs, [('ry', 'rx', 'toward l = 90° (ly)', 'toward the galactic centre, downward (ly)', 'Seen from above the galactic plane'),
                                          ('ry', 'rz', 'toward l = 90° (ly)', 'north of the galactic plane (ly)', 'Turned 90°: seen edge-on')]):
    frame(ax, 19, xl, yl)
    ia, ib = {'rx': 0, 'ry': 1, 'rz': 2}[a], {'rx': 0, 'ry': 1, 'rz': 2}[b]
    for s in core:
        dv = s['v'] - v0; c = colour(s['cls'])
        ax.plot([s[a], s[a] + dv[ia]*MY*KMS], [s[b], s[b] + dv[ib]*MY*KMS], color=c, lw=1.4, alpha=.9)
        ax.scatter([s[a]], [s[b]], s=size(s['V'])*1.6, color=c, edgecolors='none', zorder=3)
        ax.annotate(s['name'], (s[a], s[b]), xytext=(6, 4), textcoords='offset points', color='#eee', fontsize=10)
    ax.plot(0, 0, '+', ms=14, color='#9fb6d8'); ax.text(.5, -1.5, 'group centre', color='#9fb6d8', fontsize=9)
    ax.set_title(ttl, color='white', fontsize=12.5)
axs[0].invert_yaxis()   # galactic centre down, as in the other maps
fig.suptitle('The Ursa Major group nucleus in its own frame: motions relative to the group’s average, over one million years', color='white', fontsize=14)
fig.text(.5, .015, 'The seven stars of the nucleus (Merak, Phecda, Megrez, Alioth, Mizar, Alcor, 78 UMa), about 80 ly from the Sun. The lines are each star’s velocity minus the '
         'group average (1–4 km/s, close to the measuring uncertainty of about 1–2 km/s). They point in scattered directions: the group is slowly dispersing, '
         'with no sign of rotation or spiral structure.', color='#9aa', fontsize=9, ha='center', wrap=True)
fig.savefig(os.path.join(here, 'ursa_major_nucleus_own_frame.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)
print('done')
