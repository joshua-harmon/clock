# How the length of Sirius's yearly cycles has changed, 5000 BC to AD 3000 (data from make_sothic_series.py)
import csv, os, math
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker
here = os.path.dirname(os.path.abspath(__file__)); BG = '#05070d'
r = np.array([[float(x) for x in row] for row in list(csv.reader(open(os.path.join(here, 'sirius_cycles_series.csv'))))[1:]])
y = r[:, 0]
def slope(col, W=100):
    return np.array([np.polyfit(y[m] - yy, r[m, col], 1)[0] for yy in y for m in [np.abs(y - yy) <= W]])
def jan1_julian(Y):     # JD of 1 January (Julian calendar), astronomical year numbering
    yy = Y + 4800 - 1; return 1 + (153*10 + 2)//5 + 365*yy + yy//4 - 32083 - 0.5
sel = (y >= -5000) & (y <= 2994)
lines = [(1, 'Heliacal rising at Memphis (the Egyptian “Sothic” year): Sirius first seen in the dawn', '#ffd27a', 2.6),
         (2, 'Sirius on the meridian at midnight (same clock time each year)', '#7fb2ff', 2.0),
         (3, 'Sun passes Sirius (conjunction in ecliptic longitude)', '#c7a6ff', 2.0),
         (4, 'Vernal-equinox (tropical) year, for reference', '#9aa', 1.4)]
fig = plt.figure(figsize=(18, 12.5), dpi=125); fig.patch.set_facecolor(BG)
ax = fig.add_axes([0.07, 0.42, 0.90, 0.50]); bx = fig.add_axes([0.07, 0.08, 0.90, 0.25])
for a_ in (ax, bx):
    a_.set_facecolor(BG); a_.tick_params(colors='#9aa', labelsize=10); a_.grid(color='#1c2636', lw=.6)
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
TK = [1 - b for b in range(5000, 0, -500)] + list(range(500, 3001, 500))
TL = ['%d BC' % b for b in range(5000, 0, -500)] + ['AD %d' % t for t in range(500, 3001, 500)]
res = {}
for col, lab, c, lw in lines:
    s = slope(col); res[col] = s; ax.plot(y[sel], s[sel], color=c, lw=lw, label=lab)
ax.axhline(365.25, color='#ff8f6b', lw=1.2, ls='--'); ax.text(-4950, 365.25 - 0.0008, '365.25 days (Julian year)', color='#ff8f6b', fontsize=10.5)
ax.axhline(365.2425, color='#6b7d99', lw=1, ls=':'); ax.text(-4950, 365.2425 + 0.0003, '365.2425 days (Gregorian calendar)', color='#8899b3', fontsize=10)
for yr, txt in ((-2599, 'Great Pyramid\n(c. 2600 BC)'), (-1320, 'Era of Menophres\n(1321 BC)'), (139, 'Censorinus: rising\non 1 Thoth (AD 139)'), (2026, 'today')):
    ax.axvline(yr, color='#445', lw=.8, ls=':'); ax.text(yr + 25, 365.2575, txt, color='#bbb', fontsize=9.5, va='top')
def lab_at(col, yr, dy=0.0004):
    i = np.argmin(np.abs(y - yr)); return res[col][i]
ax.text(2950, lab_at(1, 2950) - 0.0009, '%.4f' % lab_at(1, 2026), color='#ffd27a', fontsize=10, ha='right')
ax.set_xlim(-5000, 3000); ax.set_ylim(365.2405, 365.2585)
ax.set_ylabel('length of the cycle (days)', color='#ccc', fontsize=11)
ax.set_xticks(TK); ax.set_xticklabels(TL); ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.4f'))
ax.legend(loc='lower right', facecolor=BG, edgecolor='#3a4a60', labelcolor='#ddd', fontsize=10)
ax.set_title('How long Sirius takes to come round again, 5000 BC to AD 3000', color='white', fontsize=14)
# lower panel: the rising date in the Julian calendar (a 365.25-day count): it barely moves
jd = r[:, 1]; doy = np.array([j - jan1_julian(int(yy)) for j, yy in zip(jd, y)]); doy = np.convolve(doy, np.ones(4)/4, mode='same')   # average over the leap-year cycle
bx.plot(y[sel], doy[sel] + 1, color='#ffd27a', lw=2)
for yr in (-2599, -1320, 139, 2026):
    i = np.argmin(np.abs(y - yr)); d = doy[i] + 1; day = d - 181
    bx.scatter([y[i]], [d], s=40, color='#ffd27a', zorder=3); bx.annotate('%s: July %d' % ('AD %d' % yr if yr > 0 else '%d BC' % (1 - yr), int(round(day))), (y[i], d), xytext=(-10, 14), ha='right', textcoords='offset points', color='#ffd27a', fontsize=9.5)
bx.set_xlim(-5000, 3000); bx.set_ylim(194, 206)
bx.set_yticks(range(195, 206, 2)); bx.set_yticklabels(['July %d' % (k - 181) for k in range(195, 206, 2)])
bx.set_xticks(TK); bx.set_xticklabels(TL)
bx.set_title('Date of Sirius’s first dawn appearance at Memphis, in the Julian calendar (365.25 days a year): only about 6 days of drift from 5000 BC to today', color='white', fontsize=12)
fig.suptitle('Sirius’s synodic cycles: when 365.25 days was (nearly) exact', color='white', fontsize=16)
fig.text(.07, .012, 'Each curve is the average cycle length over ±100 years, in days of the time (allowing for the slowing of the Earth’s rotation). Heliacal rising: Sirius 0.5° up while the Sun is 9° below the horizon '
         'at Memphis (30° N); this reproduces Censorinus’s date for AD 139. Sirius’s motion: Hipparcos (2007) and RV −5.5 km/s; precession: Vondrák et al. (2011) long-term model; Sun: Meeus. '
         'The rising year was within about a minute of 365.25 days from before 5000 BC to about 1300 BC, which is why the Egyptian 1,460-year Sothic cycle works.', color='#9aa', fontsize=9.5, wrap=True)
fig.savefig(os.path.join(here, 'sirius_synodic_cycle_over_time.png'), facecolor=BG, bbox_inches='tight'); plt.close(fig)
for col, lab, c, lw in lines: print(lab[:30], ' '.join('%d:%.4f' % (yy, res[col][np.argmin(np.abs(y - yy))]) for yy in (-4000, -2600, -1320, 139, 2026)))
