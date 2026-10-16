"""The Sothic (heliacal-rising) year at Memphis over a whole precession cycle.
Simple geometric model: Sirius fixed at its J2000 ecliptic latitude (-39.61 deg), its longitude carried round by precession
(50.29"/yr; one turn in about 25,800 years); circular solar orbit; obliquity 23.44 deg; Sirius first seen when it is 0.5 deg up
and the Sun 9 deg below the horizon at Memphis (29.85 N). Sirius's own motion, refraction details and changes in the Earth's orbit
are left out. Writes images/U4/strsWithin50ly/sirius_sothic_full_precession_cycle.png"""
import math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
EPS = math.radians(23.44); LAT = math.radians(29.85); HST = math.radians(0.5); HSUN = math.radians(-9.0)
LAM0 = 104.08; BET = math.radians(-39.61); PREC = 50.29/3600
def eq(lam, bet):
    l = math.radians(lam)
    de = math.asin(math.sin(bet)*math.cos(EPS) + math.cos(bet)*math.sin(EPS)*math.sin(l))
    ra = math.atan2(math.sin(l)*math.cos(EPS) - math.tan(bet)*math.sin(EPS), math.cos(l))
    return ra, de
def sunalt_at_rising(ls, lsun):
    ra, de = eq(ls, BET); c = (math.sin(HST) - math.sin(LAT)*math.sin(de))/(math.cos(LAT)*math.cos(de))
    if abs(c) > 1: return float('nan')
    lst = ra - math.acos(c); sra, sde = eq(lsun, 0.0)
    return math.asin(math.sin(LAT)*math.sin(sde) + math.cos(LAT)*math.cos(sde)*math.cos(lst - sra))
def offset(t):
    ls = LAM0 + PREC*t; prev = None
    for d in np.arange(-40, 140, 0.02):
        a = sunalt_at_rising(ls, ls + d)
        if a != a: return float('nan')
        if prev is not None and prev[1] > HSUN >= a: return prev[0] + (prev[1] - HSUN)/(prev[1] - a)*0.02
        prev = (d, a)
    return float('nan')
P = 360/PREC
yrs = np.arange(-24000, 28001, 25); t = yrs - 2000.0
D = np.array([offset(x) for x in t]); dec = np.array([math.degrees(eq(LAM0 + PREC*x, BET)[1]) for x in t])
date = 0.00637*t + D*365.2422/360
L = 365.25 + np.gradient(date, yrs)
L[np.isnan(D)] = np.nan
gap = yrs[np.isnan(D)]
runs = []; s0 = None
for y0, y1 in zip(yrs[:-1], yrs[1:]):
    a_, b_ = y0 in gap, y1 in gap
    if b_ and not a_: s0 = y1
    if a_ and not b_ and s0 is not None: runs.append((s0, y0)); s0 = None
print('precession period %.0f yr; no rising at Memphis:' % P, runs)
ok = ~np.isnan(L); print('cycle length range %.4f to %.4f; mean over cycle = sidereal' % (np.nanmin(L[ok & (abs(L - 365.25) < .02)]), np.nanmax(L[ok & (abs(L - 365.25) < .02)])))
def lab(y): return '%d BC' % (1 - y) if y <= 0 else 'AD %d' % y
BG = '#05070d'
fig = plt.figure(figsize=(17, 10), dpi=125); fig.patch.set_facecolor(BG)
ax = fig.add_axes([0.07, 0.42, 0.88, 0.48]); bx = fig.add_axes([0.07, 0.10, 0.88, 0.22], sharex=ax)
for a_ in (ax, bx):
    a_.set_facecolor('#0a0f1c')
    for sp in a_.spines.values(): sp.set_color('#3a4a60')
    a_.tick_params(colors='#9aa', labelsize=9)
    for g0, g1 in runs: a_.axvspan(g0, g1, color='#ff6b6b', alpha=.12, lw=0)
    a_.axvspan(-5000, 3000, color='#ffd36e', alpha=.07, lw=0)
    a_.axvline(2026, color='#ddd', lw=.8, ls=(0, (3, 3)))
ax.plot(yrs, L, color='#ffd36e', lw=2.2)
ax.axhline(365.25, color='#ff9b6b', lw=1, ls=(0, (5, 4))); ax.axhline(365.25636, color='#c9a6ff', lw=1, ls=(0, (2, 3)))
ax.axhline(365.2422, color='#9aa7b0', lw=1, ls=(0, (1, 3)))
ax.set_ylim(365.240, 365.268); ax.ticklabel_format(axis='y', useOffset=False)
ax.text(-23800, 365.2504, '365.25 days (Julian year)', color='#ff9b6b', fontsize=10)
ax.text(-23800, 365.2569, 'sidereal year 365.2564 (the average over a whole cycle)', color='#c9a6ff', fontsize=10)
ax.text(-23800, 365.2426, 'tropical year 365.2422', color='#9aa7b0', fontsize=10)
for g0, g1 in runs: ax.text((g0 + g1)/2, 365.2665, 'Sirius never rises\nat Memphis', color='#ff9b9b', fontsize=10, ha='center', va='top')
ax.text(-1000, 365.2665, 'span of the\n5000 BC – AD 3000\nchart', color='#ffd36e', fontsize=9.5, ha='center', va='top')
ax.text(2300, 365.2405, 'today', color='#ddd', fontsize=9.5, va='bottom')
ax.set_ylabel('length of the heliacal-rising year (days)', color='#ccc', fontsize=10)
ax.set_title('How long Sirius takes to come round again at Memphis, over a whole turn of precession (about %s years)' % format(int(round(P, -2)), ','), color='white', fontsize=12.5)
bx.plot(yrs, dec, color='#8fd3ff', lw=2); bx.axhline(-(90 - 29.85), color='#ff9b9b', lw=1, ls=(0, (4, 3)))
bx.text(-23800, -58.5, 'below this declination Sirius stays under the horizon at Memphis', color='#ff9b9b', fontsize=9.5)
bx.set_ylabel('Sirius’s declination (°)', color='#ccc', fontsize=10); bx.set_ylim(-66, -12)
tk = [1 - n for n in range(24000, 0, -4000)] + list(range(4000, 28001, 4000)); bx.set_xticks(tk); bx.set_xticklabels([lab(y) for y in tk])
plt.setp(ax.get_xticklabels(), visible=False); bx.set_xlim(-24000, 28000)
bx.set_title('Precession carries Sirius north and south; it returns after one full turn', color='white', fontsize=12)
fig.text(0.5, 0.955, 'Sirius’s rising year over a whole precession cycle: the arc closes after about 25,800 years', color='white', fontsize=15, ha='center')
fig.text(0.01, 0.005, 'Simple geometric model: Sirius’s latitude from the ecliptic held fixed and its longitude carried round by precession (50.29″ a year); circular solar orbit; obliquity 23.44°; '
         'Sirius first seen 0.5° up with the Sun 9° down at Memphis (29.85° N).\nIt leaves out Sirius’s own motion (about 9° over the cycle), the slow changes in the Earth’s orbit and tilt, and refraction details, so it differs from the detailed '
         '5000 BC – AD 3000 chart by up to about 0.001 day; the far ends are less certain.', color='#9aa', fontsize=9.5)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'images', 'U4', 'strsWithin50ly', 'sirius_sothic_full_precession_cycle.png')
fig.savefig(out, facecolor=BG, bbox_inches='tight'); print('saved')
