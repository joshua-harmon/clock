"""The whole galaxy as a mix of many stars, each on its own rosette path.
Test-particle model: flat rotation curve + rotating bar + rotating two-armed spiral, seen from above (north galactic pole; galaxy turns clockwise).
Writes time/u6/u6_docs/pleiades/galaxy_many_rosettes.png"""
import math, os, sys, time
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
N = int(sys.argv[1]) if len(sys.argv) > 1 else 150000
T = 1000.0; DT = 0.5                                  # million years
KMS = 1.0227e-3                                       # km/s -> kpc per million years
V0 = 230*KMS; RC = 0.3; R0 = 8.15; H = 2.6            # flat rotation curve 230 km/s; Sun 8.15 kpc; disk scale length 2.6 kpc
OB = 39*KMS; RB = 5.0; ALPHA = 0.012                  # bar: pattern speed 39 km/s/kpc, radius 5 kpc, strength (Dehnen 2000 form)
OS = 22*KMS; M = 2; PITCH = math.radians(12); EPS = 0.012   # spiral: pattern speed 22 km/s/kpc, 2 arms, pitch 12 deg
AB = ALPHA*V0**2/3*(R0/RB)**3
KS = M/math.tan(PITCH)
BAR_END = math.radians(90 - 27)                      # near end of the bar 27 deg from the Sun-centre line, toward l = 0..90
PHB0 = BAR_END + OB*T
CHI0 = math.pi - M*(math.pi/2 + OS*T)                # Sun (at angle 90 deg) between arms at the end
rng = np.random.default_rng(7)
RMAX = 20.0
half = N//2
Ra = rng.gamma(2.0, H, size=3*N); Ra = Ra[(Ra > 0.4) & (Ra < RMAX)][:half]          # half drawn from the exponential disk
Rb = RMAX*np.sqrt(rng.uniform((0.4/RMAX)**2, 1, N - half))                          # half spread evenly in area (better coverage far out)
R = np.concatenate([Ra, Rb])
_pa = R*np.exp(-R/H)/(H**2); _pb = 2*R/RMAX**2                                       # the two sampling densities (per unit R)
WEXP = R*np.exp(-R/H)/(0.5*_pa + 0.5*_pb)                                          # weights that restore the exponential disk
ph = rng.uniform(0, 2*math.pi, N)
vc = V0*R/np.sqrt(R**2 + RC**2)
sR = np.minimum(40*np.exp(-(R - R0)/(2*H)), 90)*KMS; sP = sR/math.sqrt(2)
vR = rng.normal(0, 1, N)*sR; vP = -(vc - (sR**2/vc)*0.5 + rng.normal(0, 1, N)*sP)   # clockwise
x, y = R*np.cos(ph), R*np.sin(ph)
vx = vR*np.cos(ph) - vP*np.sin(ph); vy = vR*np.sin(ph) + vP*np.cos(ph)
def acc(x, y, t):
    R2 = x*x + y*y; R = np.sqrt(R2); ph = np.arctan2(y, x)
    aR = -V0**2*R/(R2 + RC**2); aP = np.zeros_like(R)
    g = min(1.0, t/300.0); g = g*g*(3 - 2*g)          # grow bar and arms smoothly over 300 Myr
    c2 = np.cos(2*(ph + OB*t - PHB0)); s2 = np.sin(2*(ph + OB*t - PHB0))
    out = R >= RB
    f = np.where(out, -(RB/R)**3, (R/RB)**3 - 2); fp = np.where(out, 3*RB**3/R**4, 3*R**2/RB**3)
    aR += -g*AB*fp*c2; aP += 2*g*AB*f*s2/R
    chi = M*(ph + OS*t) - KS*np.log(R/R0) + CHI0
    sc = np.sin(chi)
    aR += -g*EPS*V0**2*sc*(-KS/R); aP += -g*EPS*V0**2*sc*M/R
    return aR*np.cos(ph) - aP*np.sin(ph), aR*np.sin(ph) + aP*np.cos(ph)
NPZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'galaxy_test_particles_end.npz')
if len(sys.argv) > 2 and sys.argv[2] == 'plot':
    d = np.load(NPZ); x, y, W = d['x'], d['y'], d['w']; N = len(x); nst = 0
t0 = time.time(); t = 0.0
ax_, ay_ = acc(x, y, t)
nst = int(T/DT) if not (len(sys.argv) > 2 and sys.argv[2] == 'plot') else 0
for k in range(nst):
    vx += 0.5*DT*ax_; vy += 0.5*DT*ay_
    x += DT*vx; y += DT*vy; t += DT
    ax_, ay_ = acc(x, y, t)
    vx += 0.5*DT*ax_; vy += 0.5*DT*ay_
print('integrated %d stars for %.0f Myr in %.0f s' % (N, T, time.time() - t0))
if nst: W = WEXP; np.savez_compressed(NPZ, x=x, y=y, w=W)

# ---------------- picture ----------------
LY = 3261.56; BG = '#05070d'
X, Y = x*LY/1000, y*LY/1000                      # thousands of light years
L = 61; NB = 610
Hd, xe, ye = np.histogram2d(X, Y, bins=NB, range=[[-L, L], [-L, L]], weights=W)
Hn, _, _ = np.histogram2d(X, Y, bins=NB, range=[[-L, L], [-L, L]]); Hn = Hn.T
Hd = Hd.T
xc = 0.5*(xe[1:] + xe[:-1]); XX, YY = np.meshgrid(xc, xc); RR = np.hypot(XX, YY)
def gaussian_filter(a, sig):
    n = a.shape[0]; k = np.fft.fftfreq(n)
    g = np.exp(-2*(math.pi*sig)**2*(k[:, None]**2 + k[None, :]**2))
    return np.real(np.fft.ifft2(np.fft.fft2(a)*g))
Hs = gaussian_filter(Hd, 1.8); Hc = gaussian_filter(Hn, 5.0)
rb = np.linspace(0, L*1.42, 140); idx = np.digitize(RR, rb)
mean = np.zeros_like(Hs)
for i in np.unique(idx):
    m_ = idx == i; mean[m_] = Hc[m_].mean()
con = np.where(mean > 0, Hc/np.maximum(mean, 1e-9), np.nan)
con[RR > 41] = np.nan; con[RR < 3] = np.nan
fig = plt.figure(figsize=(21, 9.8), dpi=125); fig.patch.set_facecolor(BG)
# Lindblad resonances for a flat rotation curve: Omega = V/R, kappa = sqrt(2) Omega
VKMS = 230.0; F = 1 - 1/math.sqrt(2); G = 1 + 1/math.sqrt(2)
RES = []
for pat, om, col in (('bar', 39.0, '#7fe3ff'), ('spiral arms', 22.0, '#9dff8a')):
    for nm, fac, ls in (('inner Lindblad resonance', F, (0, (1.5, 2.5))), ('corotation', 1.0, '-'), ('outer Lindblad resonance', G, (0, (7, 4)))):
        rk = VKMS*fac/om; RES.append((pat, nm, rk*LY/1000, col, ls))
        print('%s %s: %.2f kpc = %.1f thousand ly' % (pat, nm, rk, rk*LY/1000))
sun = (0, R0*LY/1000)
for k in range(2):
    ax = fig.add_axes([0.015 + 0.40*k, 0.12, 0.385, 0.76]); ax.set_facecolor('black')
    if k == 0:
        img = np.log10(Hs + Hs.max()*1e-4); ax.imshow(img, extent=[-L, L, -L, L], origin='lower', cmap='afmhot', vmin=np.percentile(img, 35), vmax=np.percentile(img, 99.7))
        ax.set_title('All the stars together: how the galaxy looks from above', color='white', fontsize=12.5)
    else:
        ax.imshow(con, extent=[-L, L, -L, L], origin='lower', cmap='magma', vmin=0.6, vmax=1.5)
        ax.set_title('Same stars, brightness compared with the average at each distance from the centre', color='white', fontsize=12.5)
    th = np.linspace(0, 2*math.pi, 400); ax.plot(R0*LY/1000*np.cos(th), R0*LY/1000*np.sin(th), color='#9fd0ff', lw=.7, ls=(0, (4, 4)), alpha=.6)
    ax.plot(*sun, marker='o', ms=8, color='#ffe24a', mec='black', mew=.8, zorder=5)
    ax.annotate('Sun', sun, xytext=(10, 8), textcoords='offset points', color='#ffe24a', fontsize=11, weight='bold')
    for pat, nm, rl, col, ls in RES:
        ax.plot(rl*np.cos(th), rl*np.sin(th), color=col, lw=1.4 if nm == 'corotation' else 1.2, ls=ls, alpha=.95, zorder=4)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect('equal')
    ax.set_xticks(range(-50, 51, 25)); ax.set_yticks(range(-50, 51, 25))
    ax.tick_params(colors='#9aa', labelsize=9)
    for sp in ax.spines.values(): sp.set_color('#3a4a60')
    ax.set_xlabel('thousands of light years', color='#ccc', fontsize=10)
fig.text(0.41, 0.955, 'The galaxy as %s stars, each on its own rosette path' % format(N, ','), color='white', fontsize=15, ha='center')
fig.text(0.41, 0.918, 'No star stays in an arm. The arms and the bar are places where the rosette paths crowd together, like a traffic jam that stays put while the cars move through it.',
         color='#e8e8e8', fontsize=11, ha='center')
yk = 0.86
fig.text(0.815, yk, 'Lindblad resonances', color='white', fontsize=13, weight='bold', va='top'); yk -= 0.04
for pat, col in (('bar', '#7fe3ff'), ('spiral arms', '#9dff8a')):
    fig.text(0.815, yk, 'the %s (turning %s km/s/kpc)' % (pat, '39' if pat == 'bar' else '22'), color=col, fontsize=11, weight='bold', va='top'); yk -= 0.035
    for p_, nm, rl, c_, ls in [r for r in RES if r[0] == pat]:
        kx = fig.add_axes([0.815, yk - 0.022, 0.03, 0.02]); kx.axis('off'); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
        kx.plot([0, 1], [0.5, 0.5], color=c_, lw=1.6, ls=ls)
        fig.text(0.85, yk, '%s: %.1f thousand ly%s' % (nm, rl, '  (off the picture)' if rl > L else ''), color=c_, fontsize=10, va='top'); yk -= 0.032
    yk -= 0.02
fig.text(0.815, yk, 'Inner and outer Lindblad resonance:\na star’s in-and-out loop keeps exact\ntime with the turning pattern, so its\nrosette closes as the pattern sees it.\nCorotation: the star keeps pace with\nthe pattern.\n\nThe arms live between the spiral’s\ninner and outer resonances. The Sun\n(26.6 thousand ly) lies between the\nbar’s corotation and its outer\nresonance, about 6,000 ly inside\nthe latter.', color='#ccc', fontsize=10, va='top', linespacing=1.3)
fig.text(0.01, 0.005, 'Test-particle model run for 1 billion years: flat rotation curve (230 km/s), a rotating bar (pattern speed 39 km/s/kpc, half-length 5 kpc; Dehnen 2000 form) and a rotating two-armed spiral '
         '(22 km/s/kpc, pitch 12°), seen from above with the galaxy turning clockwise.\nThe bar and spiral strengths are set by hand, so this shows how arms arise, not the Milky Way’s exact arms. '
         'Real pictures of the arms also owe much to gas, dust and young bright stars, which are not in this model. Faint dashed circle = the Sun’s distance from the centre (26,600 ly).\nResonance radii assume a flat rotation curve (230 km/s). Measured bar speeds range from about 35 to 55 km/s/kpc, which moves the bar’s outer resonance from beyond the Sun to just inside it.',
         color='#9aa', fontsize=9.5)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'time', 'u6', 'u6_docs', 'pleiades', 'galaxy_many_rosettes.png'), facecolor=BG, bbox_inches='tight'); print('saved')
