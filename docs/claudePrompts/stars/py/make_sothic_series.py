import sys, csv, math, warnings; warnings.filterwarnings('ignore')
import numpy as np, erfa
import sothic_model as sm
def rising_fast(year, guess):
    jd0 = guess - 4; prev = None
    for k in range(0, 12):
        h, t = sm.sun_alt_at_rising(jd0 + k)
        if prev is not None and prev[0] > sm.H_SUN >= h:
            f = (prev[0] - sm.H_SUN)/(prev[0] - h); return prev[1] + f*(t - prev[1])
        prev = (h, t)
    return sm.heliacal_rising(year)
years = np.arange(-5200, 3201, 2); rows = []; g = None
for y in years:
    hr = sm.heliacal_rising(y) if g is None else rising_fast(y, g + 2*365.25)
    g = hr
    rows.append((y, hr, sm.midnight_transit(y), sm.conjunction(y), sm.equinox(y)))
with open('sirius_cycles_series.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['year', 'heliacal_rising_JD_UT', 'midnight_transit_JD_UT', 'conjunction_JD_UT', 'vernal_equinox_JD_UT']); w.writerows(rows)
print('done', len(rows))
