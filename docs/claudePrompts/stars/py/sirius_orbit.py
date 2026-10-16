# Sirius A's reflex orbit around the Sirius AB barycentre (Bond et al. 2017 elements)
import math, numpy as np
P, A, E, I, NODE, OM, T0 = 50.1284, 7.4957, 0.59142, 136.336, 45.400, 149.161, 1994.5715
MA, MB = 2.063, 1.018; F = MB/(MA + MB)
def rel(t):
    """B relative to A, arcsec (east, north)"""
    M = 2*math.pi*(t - T0)/P; Ea = M
    for _ in range(50): Ea = Ea - (Ea - E*math.sin(Ea) - M)/(1 - E*math.cos(Ea))
    nu = 2*math.atan2(math.sqrt(1+E)*math.sin(Ea/2), math.sqrt(1-E)*math.cos(Ea/2))
    r = A*(1 - E*math.cos(Ea)); u = math.radians(OM) + nu; O, i = math.radians(NODE), math.radians(I)
    north = r*(math.cos(O)*math.cos(u) - math.sin(O)*math.sin(u)*math.cos(i))
    east = r*(math.sin(O)*math.cos(u) + math.cos(O)*math.sin(u)*math.cos(i))
    return np.array([east, north])
def a_off(t): return -F*rel(t)          # Sirius A from the barycentre, arcsec (east, north)
def a_vel(t, h=1e-3): return (a_off(t+h) - a_off(t-h))/(2*h)   # arcsec/yr
