# Sirius over the last 100,000 years: its track across the sky (from Hipparcos proper motion and radial velocity,
# straight-line motion in 3-D) and its changing distance and brightness.
import csv, math, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
K = 4.740470446; KMS_LY = 1.0227e-6*3.26156          # km/s -> ly per year
RA0, DE0, PLX, PMRA, PMDE, RV, V0 = 101.2872, -16.7161, 379.21, -546.01, -1223.07, -5.5, -1.46   # Hipparcos 2007; RV Pulkovo

def unit(ra, de):
    a, d = math.radians(ra), math.radians(de); return np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
r0 = unit(RA0, DE0); a, d = math.radians(RA0), math.radians(DE0)
p = np.array([-math.sin(a), math.cos(a), 0]); q = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
D0 = 3261.56/PLX
vel = ((K*PMRA/PLX)*p + (K*PMDE/PLX)*q + RV*r0)*KMS_LY          # ly per year, equatorial
years = np.arange(-100000, 1, 1000)
pos = D0*r0[None, :] + years[:, None]*vel[None, :]
dist = np.linalg.norm(pos, axis=1)
ra = np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) % 360; de = np.degrees(np.arcsin(pos[:, 2]/dist))
mag = V0 + 5*np.log10(dist/D0)
print('100,000 years ago: RA %.1f Dec %.1f, %.1f ly, magnitude %.2f; moved %.1f deg' % (ra[0], de[0], dist[0], mag[0],
      math.degrees(math.acos(np.dot(pos[0]/dist[0], r0)))))

# background: naked-eye stars from the 250 ly list in the same part of the sky (their own motions are much smaller, shown fixed)
bg = []
for r in csv.DictReader(open(os.path.join(here, 'naked_eye_stars_within_250ly.csv'), encoding='utf-8')):
    h, m = r['RA (h m)'].split(); sgn = -1 if r['Dec (d m)'][0] == '-' else 1; dd, dm = r['Dec (d m)'][1:].split()
    bg.append((15*(int(h) + float(m)/60), sgn*(int(dd) + int(dm)/60), float(r['V mag']), r['Proper name'], r['Name']))
LAB = {'Procyon', 'Canopus', 'Aldebaran', 'Castor', 'Pollux', 'Alphard', 'Regulus', 'Capella'}
FAR = [('Betelgeuse', 88.79, 7.41, 0.45), ('Rigel', 78.63, -8.20, 0.13), ('Bellatrix', 81.28, 6.35, 1.64), ('Adhara', 104.66, -28.97, 1.50),
       ('Wezen', 107.10, -26.39, 1.83), ('Mirzam', 95.67, -17.96, 1.98), ('Alnilam', 84.05, -1.20, 1.69), ('Saiph', 86.94, -9.67, 2.07)]

fig = plt.figure(figsize=(18, 9), dpi=130); fig.patch.set_facecolor('#05070d')
ax = fig.add_axes([0.04, 0.1, 0.56, 0.8]); ax.set_facecolor('#05070d')
for s in ax.spines.values(): s.set_color('#3a4a60')
ax.tick_params(colors='#9aa'); ax.grid(color='#1c2636', lw=.6)
def sz(v): return max(3, (6.6 - v)**2*4)
for ra_, de_, v, pn, nm in bg:
    if 60 < ra_ < 150 and -45 < de_ < 25:
        ax.scatter([ra_], [de_], s=sz(v), color='#c9d4e6', edgecolors='none', alpha=.8)
        if pn in LAB: ax.annotate(pn, (ra_, de_), xytext=(5, 3), textcoords='offset points', color='#ccc', fontsize=8.5)
for nm, ra_, de_, v in FAR:
    ax.scatter([ra_], [de_], s=sz(v), color='#c9d4e6', edgecolors='none'); ax.annotate(nm, (ra_, de_), xytext=(5, 3), textcoords='offset points', color='#ccc', fontsize=8.5)
ax.plot(ra, de, color='#7fb2ff', lw=1.6, zorder=3)
for k in range(0, len(years), 10):
    ax.scatter([ra[k]], [de[k]], s=sz(mag[k])*0.9, color='#e6ecff', edgecolors='#7fb2ff', lw=.8, zorder=4)
    ax.annotate('%s yr ago' % format(-years[k], ',') if years[k] else 'today', (ra[k], de[k]), xytext=(7, -10), textcoords='offset points', color='#9fc4ff', fontsize=8.5)
ax.scatter([RA0], [DE0], s=sz(V0), color='#ffffff', zorder=5); ax.annotate('Sirius', (RA0, DE0), xytext=(8, 6), textcoords='offset points', color='white', fontsize=11, weight='bold')
ax.set_xlim(150, 60); ax.set_ylim(-45, 25); ax.set_aspect('equal')
ax.set_xlabel('right ascension (degrees; east to the left)', color='#ccc'); ax.set_ylabel('declination (degrees)', color='#ccc')
ax.set_title('Sirius’s path across the sky over the last 100,000 years', color='white', fontsize=13)

ax2 = fig.add_axes([0.67, 0.55, 0.31, 0.33]); ax3 = fig.add_axes([0.67, 0.1, 0.31, 0.33])
for a_ in (ax2, ax3):
    a_.set_facecolor('#05070d'); a_.tick_params(colors='#9aa'); a_.grid(color='#1c2636', lw=.6)
    for s in a_.spines.values(): s.set_color('#3a4a60')
ax2.plot(years/1000, dist, color='#7fb2ff', lw=1.8); ax2.set_ylabel('distance from the Sun (ly)', color='#ccc'); ax2.set_title('Distance', color='white', fontsize=11)
ax3.plot(years/1000, mag, color='#ffd27a', lw=1.8); ax3.invert_yaxis(); ax3.set_ylabel('apparent magnitude (brighter is up)', color='#ccc')
ax3.set_xlabel('thousands of years from now', color='#ccc'); ax3.set_title('Brightness seen from Earth', color='white', fontsize=11)
fig.suptitle('Sirius over the last 100,000 years: moved %.0f° across the sky, %.1f ly away then, %.1f ly now' %
             (math.degrees(math.acos(np.dot(pos[0]/dist[0], r0))), dist[0], D0), color='white', fontsize=14)
fig.text(.04, .02, 'Straight-line motion from today’s Hipparcos position, parallax and proper motion (1.34″ a year) and a radial velocity of −5.5 km/s (Pulkovo). '
         'Marks every 10,000 years. The other stars are shown at today’s places; most move far less, but Procyon also travels fast (1.26″ a year, about 35° in 100,000 years). Precession of the equinoxes is not applied: these are fixed J2000 star coordinates.',
         color='#9aa', fontsize=9, wrap=True)
fig.savefig(os.path.join(here, 'sirius_motion_100000_years.png'), facecolor=fig.get_facecolor(), bbox_inches='tight'); plt.close(fig)
print('done')
