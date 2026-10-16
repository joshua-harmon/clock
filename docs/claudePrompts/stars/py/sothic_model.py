# Sirius cycles over 6,000 years: heliacal rising at Memphis (the Egyptian "Sothic" year), midnight meridian transit,
# and conjunction with the Sun. Sirius: Hipparcos position + 3-D space motion; precession: Vondrak et al. 2011 long-term
# model (erfa.ltpb); Sun: Meeus low-precision theory (mean equinox of date); Earth-rotation: Morrison & Stephenson dT.
import math, numpy as np, erfa
K = 4.740470446; KMS_LY = 1.0227e-6*3.26156
RA0, DE0, PLX, PMRA, PMDE, RV = 101.2872, -16.7161, 379.21, -546.01, -1223.07, -5.5
a, d = math.radians(RA0), math.radians(DE0)
r0 = np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
p = np.array([-math.sin(a), math.cos(a), 0]); q = np.array([-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d)])
D0 = 3261.56/PLX; VEL = ((K*PMRA/PLX)*p + (K*PMDE/PLX)*q + RV*r0)*KMS_LY
def jan1(year): return 2451544.5 + (year - 2000)*365.25     # Julian-calendar 1 January (approx.)
def jd_to_epj(jd): return 2000.0 + (jd - 2451545.0)/365.25
def deltaT(jd):  # seconds
    u = (jd_to_epj(jd) - 1820)/100; return -20 + 32*u*u
def sirius_radec(jd_tt):
    ep = jd_to_epj(jd_tt); v = D0*r0 + (ep - 1991.25)*VEL; v = erfa.ltpb(ep) @ (v/np.linalg.norm(v))
    return math.atan2(v[1], v[0]) % (2*math.pi), math.asin(v[2])
def obliquity(jd_tt): return erfa.obl06(jd_tt, 0.0)
def sun_lon(jd_tt):  # apparent-ish ecliptic longitude, mean equinox of date, radians
    T = (jd_tt - 2451545.0)/36525
    L0 = 280.46646 + 36000.76983*T + 0.0003032*T*T; M = math.radians(357.52911 + 35999.05029*T - 0.0001537*T*T)
    C = (1.914602 - 0.004817*T - 0.000014*T*T)*math.sin(M) + (0.019993 - 0.000101*T)*math.sin(2*M) + 0.000289*math.sin(3*M)
    return math.radians((L0 + C - 0.00569) % 360)
def sun_radec(jd_tt):
    l, e = sun_lon(jd_tt), obliquity(jd_tt)
    return math.atan2(math.cos(e)*math.sin(l), math.cos(l)) % (2*math.pi), math.asin(math.sin(e)*math.sin(l))
def gmst(jd_ut): return math.radians((280.46061837 + 360.98564736629*(jd_ut - 2451545.0)) % 360)
def ecl_lon(ra, de, e): return math.atan2(math.sin(ra)*math.cos(e) + math.tan(de)*math.sin(e), math.cos(ra)) % (2*math.pi)

LAT, LON = math.radians(29.85), math.radians(31.25)       # Memphis
H_STAR, H_SUN = math.radians(0.5), math.radians(-9.0)     # Sirius just up (about 1° apparent), Sun 9° below: Sirius first seen
def sun_alt_at_rising(jd_day):          # jd_day: UT noon-based day number; returns Sun's altitude when Sirius rises that morning
    jd = jd_day
    for _ in range(3):
        ra, de = sirius_radec(jd + deltaT(jd)/86400)
        H0 = math.acos((math.sin(H_STAR) - math.sin(LAT)*math.sin(de))/(math.cos(LAT)*math.cos(de)))
        lst_target = ra - H0
        dl = (lst_target - (gmst(jd) + LON) + math.pi) % (2*math.pi) - math.pi
        jd = jd + dl/(2*math.pi*1.00273790935)
    sra, sde = sun_radec(jd + deltaT(jd)/86400); H = gmst(jd) + LON - sra
    return math.asin(math.sin(LAT)*math.sin(sde) + math.cos(LAT)*math.cos(sde)*math.cos(H)), jd
def heliacal_rising(year):              # continuous (fractional) UT JD of first visibility in the given (astronomical) year
    jd0 = jan1(year) + 0.25   # dawn-ish start
    guess = jd0 + 120 + (-(year - 139)*0.0) ; prev = None
    for k in range(90, 330):
        h, t = sun_alt_at_rising(jd0 + k)
        if prev is not None and prev[0] > H_SUN >= h:
            f = (prev[0] - H_SUN)/(prev[0] - h); return prev[1] + f*(t - prev[1])
        prev = (h, t)
    return float('nan')
def solve(f, target, jd):               # Newton-ish root of angle f(jd) = target, TT
    for _ in range(30):
        err = (f(jd) - target + math.pi) % (2*math.pi) - math.pi
        if abs(err) < 1e-10: break
        jd -= err/(2*math.pi/365.2422)
    return jd
def midnight_transit(year):             # Sun's RA = Sirius's RA + 180°: Sirius on the meridian at (apparent) midnight; returns UT JD
    jd = jan1(year)
    jd = solve(lambda t: sun_radec(t)[0] - sirius_radec(t)[0] - math.pi, 0, jd)
    return jd - deltaT(jd)/86400
def conjunction(year):                  # Sun's ecliptic longitude = Sirius's
    jd = jan1(year) + 181
    def f(t):
        ra, de = sirius_radec(t); return sun_lon(t) - ecl_lon(ra, de, obliquity(t))
    jd = solve(f, 0, jd); return jd - deltaT(jd)/86400
def equinox(year):
    jd = jan1(year) + 79
    jd = solve(sun_lon, 0, jd); return jd - deltaT(jd)/86400
