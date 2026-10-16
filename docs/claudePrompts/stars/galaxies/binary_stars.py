# Shared data and orbit model for the two nearby white-dwarf binaries, Sirius and Procyon.
# Barycentre: Hipparcos (van Leeuwen 2007) position, parallax and proper motion + radial velocity, straight-line motion.
# Orbits: Sirius, Bond et al. 2017; Procyon, Bond et al. 2015 (via the Wikipedia star pages).
import math, numpy as np
T = np.array([[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]])
K = 4.740470446; KMS_LY = 1.0227e-6*3.26156; AU_LY = 63241.077; EP0 = 1991.25
STARS = {
 'Sirius':  dict(ra=101.2872, de=-16.7161, plx=379.21, pmra=-546.01, pmde=-1223.07, rv=-5.5,
                 P=50.1284, a=7.4957, e=0.59142, i=136.336, node=45.400, om=149.161, T0=1994.5715, mA=2.063, mB=1.018,
                 hist=[(1718, 'Halley: first proper motions'), (1844, 'Bessel infers a companion'), (1862, 'Alvan Clark sees Sirius B')],
                 close=(1830, 1880), colour='#7fb2ff', B='Sirius B'),
 'Procyon': dict(ra=114.8255, de=5.2250, plx=284.56, pmra=-714.59, pmde=-1036.80, rv=-3.2,
                 P=40.840, a=4.3075, e=0.39785, i=31.408, node=100.683, om=89.23, T0=1968.076, mA=1.478, mB=0.592,
                 hist=[(1718, 'Halley: first proper motions'), (1844, 'Bessel infers a companion'), (1862, 'Auwers computes the orbit'), (1896, 'Schaeberle sees Procyon B')],
                 close=(1856, 1897), colour='#ffe08a', B='Procyon B'),
}
class Binary:
    def __init__(s, name):
        s.name = name; s.__dict__.update(STARS[name])
        a, d = math.radians(s.ra), math.radians(s.de)
        s.r0 = np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
        s.p = np.array([-math.sin(a), math.cos(a), 0]); s.q = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
        s.D0 = 3261.56/s.plx; s.vel = ((K*s.pmra/s.plx)*s.p + (K*s.pmde/s.plx)*s.q + s.rv*s.r0)*KMS_LY
        s.F = s.mB/(s.mA + s.mB); s.au = 1000/s.plx            # AU per arcsec
    def rel(s, t):       # B relative to A: (east", north", away")
        M = (2*math.pi*(t - s.T0)/s.P) % (2*math.pi); Ea = M + s.e*math.sin(M)
        for _ in range(60): Ea -= (Ea - s.e*math.sin(Ea) - M)/(1 - s.e*math.cos(Ea))
        nu = 2*math.atan2(math.sqrt(1+s.e)*math.sin(Ea/2), math.sqrt(1-s.e)*math.cos(Ea/2))
        r = s.a*(1 - s.e*math.cos(Ea)); u = math.radians(s.om) + nu; O, i = math.radians(s.node), math.radians(s.i)
        return np.array([r*(math.sin(O)*math.cos(u) + math.cos(O)*math.sin(u)*math.cos(i)),
                         r*(math.cos(O)*math.cos(u) - math.sin(O)*math.sin(u)*math.cos(i)), r*math.sin(u)*math.sin(i)])
    def a_off(s, t): return -s.F*s.rel(t)[:2]                 # A from barycentre on the sky, arcsec
    def sky_bary(s, y):                                       # arcsec east/north of a fixed reference
        v = s.D0*s.r0 + (y - EP0)*s.vel; u = v/np.linalg.norm(v); return 206264.806*np.array([np.dot(u, s.p), np.dot(u, s.q)])
    def sky_A(s, y): return s.sky_bary(y) + s.a_off(y)
    # 3-D, galactic, AU from the Sun
    def rel3(s, t): e, n, w = s.rel(t)*s.au; return T @ (e*s.p + n*s.q + w*s.r0)
    def bary3(s, t): return T @ ((s.D0*s.r0 + (t - EP0)*s.vel)*AU_LY)
    def A3(s, t): return s.bary3(t) - s.F*s.rel3(t)
    def B3(s, t): return s.bary3(t) + (1 - s.F)*s.rel3(t)
    def speeds(s):
        g = T @ s.vel; return np.linalg.norm(s.vel)/KMS_LY, math.hypot(g[0], g[1])/KMS_LY, g[2]/KMS_LY
def top(v): return np.array([v[1], -v[0]])
