# Side-by-side comparison of the two nearby white-dwarf binaries, Sirius and Procyon, at matching scales.
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from binary_stars import Binary, top, T, KMS_LY, AU_LY

here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
S, Pr = Binary('Sirius'), Binary('Procyon')
def style(ax, title, fs=12):
    ax.set_facecolor(BG); ax.tick_params(colors='#9aa', labelsize=9); ax.grid(color='#1c2636', lw=.6)
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.set_title(title, color='white', fontsize=fs)

fig = plt.figure(figsize=(21, 13.5), dpi=120); fig.patch.set_facecolor(BG)
axs = [fig.add_axes([0.03 + k*0.21, 0.52, 0.19, 0.36]) for k in range(2)]
axg = [fig.add_axes([0.03 + k*0.21, 0.07, 0.19, 0.36]) for k in range(2)]
axw = fig.add_axes([0.47, 0.52, 0.51, 0.36]); axt = fig.add_axes([0.47, 0.07, 0.51, 0.36]); axt.axis('off')

for ax, b in zip(axs, (S, Pr)):        # on the sky, arcseconds, same scale
    style(ax, '%s on the sky (arcseconds)' % b.name)
    ty = np.arange(0, b.P, 0.05) + b.T0
    oa = np.array([b.a_off(t) for t in ty]); ob = np.array([b.a_off(t) + b.rel(t)[:2] for t in ty])
    ax.plot(ob[:, 0], ob[:, 1], color='#bbb', lw=1, label=b.B); ax.plot(oa[:, 0], oa[:, 1], color=b.colour, lw=2, label='%s A' % b.name)
    ax.scatter([0], [0], marker='+', s=70, color='#ffe24a', zorder=4)
    va, vb = b.a_off(2026), b.a_off(2026) + b.rel(2026)[:2]
    ax.plot([va[0], vb[0]], [va[1], vb[1]], color='white', lw=.7, ls=':'); ax.scatter([va[0], vb[0]], [va[1], vb[1]], s=[40, 16], color='white', zorder=5)
    ax.annotate('2026', vb, xytext=(5, 4), textcoords='offset points', color='white', fontsize=9)
    ax.set_xlim(8.5, -8.5); ax.set_ylim(-8.5, 8.5); ax.set_aspect('equal')
    ax.set_xlabel('east ←   arcsec   → west', color='#ccc'); ax.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=8.5)
for ax, b in zip(axg, (S, Pr)):        # from the north galactic pole, AU, same scale
    style(ax, '%s from the north galactic pole (AU)' % b.name)
    to = np.arange(0, b.P, 0.02) + b.T0
    oa = np.array([top(b.A3(t) - b.bary3(t)) for t in to]); ob = np.array([top(b.B3(t) - b.bary3(t)) for t in to])
    ax.plot(ob[:, 0], ob[:, 1], color='#bbb', lw=1); ax.plot(oa[:, 0], oa[:, 1], color=b.colour, lw=2); ax.scatter([0], [0], marker='+', s=70, color='#ffe24a', zorder=4)
    mv = top(T @ b.vel); mv /= np.linalg.norm(mv)
    ax.annotate('', xy=(-17 + 6*mv[0], -17 + 6*mv[1]), xytext=(-17, -17), arrowprops=dict(arrowstyle='-|>', color=b.colour, lw=1.6))
    ax.text(-17, -12.5 if mv[1] < 0 else -20.5, 'travel', color=b.colour, fontsize=8.5, ha='center')
    ax.set_xlim(-24, 24); ax.set_ylim(-24, 24); ax.set_aspect('equal'); ax.set_xlabel('AU toward l = 90° →   (galactic centre down)', color='#ccc')

# along-track wobble over time, both stars
style(axw, 'How far each bright star runs ahead of (+) or behind (−) its steady motion, 1700 to 2026 (arcseconds)')
yrs = np.arange(1700, 2026.01, 0.1)
for b in (S, Pr):
    d = np.array([b.pmra, b.pmde]); d /= np.linalg.norm(d)
    axw.plot(yrs, [b.a_off(y) @ d for y in yrs], color=b.colour, lw=1.8, label='%s A (orbit %.1f years)' % (b.name, b.P))
axw.axhline(0, color='#556', lw=.8)
for y, nm, ty, ha in ((1844, 'Bessel infers both companions', 4.5, 'right'), (1862, 'Clark sees Sirius B', 4.5, 'left'), (1896, 'Schaeberle sees Procyon B', 3.85, 'left')):
    axw.axvline(y, color='#ffb36b', lw=.7, ls=':'); axw.text(y + (1.5 if ha == 'left' else -1.5), ty, '%d: %s' % (y, nm), color='#ffb36b', fontsize=9, va='top', ha=ha)
axw.set_xlim(1700, 2026); axw.set_ylim(-2.2, 4.7); axw.set_xlabel('year', color='#ccc'); axw.set_ylabel('arcseconds', color='#ccc')
axw.legend(loc='lower left', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=9.5, ncol=2)

# comparison table
def row(b):
    sp, vin, vz = b.speeds(); wob = max(math.hypot(*b.a_off(b.T0 + k*b.P/400)) for k in range(400))
    trk = math.hypot(*(b.sky_A(1526) - b.sky_A(2026)))/60
    return [b.name, '%.2f ly' % b.D0, '%.2f + %.2f' % (b.mA, b.mB), '%.2f yr' % b.P, '%.2f' % b.e, '%.1f – %.1f AU' % (b.a*b.au*(1 - b.e), b.a*b.au*(1 + b.e)),
            '%.0f° from face-on' % (b.i if b.i < 90 else 180 - b.i), '%.3f″/yr' % (math.hypot(b.pmra, b.pmde)/1000), '%.1f′' % trk, '%.1f″' % wob, '%.1f km/s' % sp, '%.0f AU' % (vin*b.P*KMS_LY*AU_LY)]
hdr = ['', 'distance', 'masses (suns)', 'orbit', 'eccentricity', 'closest – farthest', 'tilt to sky', 'proper motion', 'moved in 500 yr', 'A’s wobble (max)', 'speed vs Sun', 'travel per orbit\n(along the plane)']
cells = list(zip(hdr, row(S), row(Pr)))
tb = axt.table(cellText=[c[1:] for c in cells], rowLabels=[c[0] for c in cells], colLabels=None, loc='upper center', cellLoc='center', colWidths=[.28, .28])
tb.auto_set_font_size(False); tb.set_fontsize(11); tb.scale(1, 1.75)
for (r, c), cell in tb.get_celld().items():
    cell.set_edgecolor('#2a3548'); cell.set_facecolor('#0b1220' if r % 2 else BG)
    cell.get_text().set_color(S.colour if (r == 0 and c == 0) else Pr.colour if (r == 0 and c == 1) else '#e6e6e6' if c >= 0 else '#9fb0c8')
    if r == 0: cell.get_text().set_weight('bold')
axt.set_title('Side by side', color='white', fontsize=12)

fig.suptitle('Sirius and Procyon compared: two bright neighbours, each pulled around by a white-dwarf companion', color='white', fontsize=15.5)
fig.text(.03, .015, 'Left panels share one scale (sky: arcseconds; space: AU) so the sizes compare directly. Sirius B is the heavier companion on a wider, more eccentric orbit, so Sirius A wobbles about 2½ times as much as Procyon A. '
         'Both pairs are heading south of the galactic plane. From Earth both orbits are seen fairly open; from the galactic pole Sirius’s is almost edge-on while Procyon’s is open. '
         'Hipparcos (2007) motions; orbits from Bond et al. 2015 (Procyon) and 2017 (Sirius).', color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, 'sirius_vs_procyon.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
print('done')
