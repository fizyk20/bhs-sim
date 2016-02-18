#!/usr/bin/env python3
import math

def delta(r):
    return 4*r**5*(r**3-27*r+54)

def cos_thresh(r):
    a = 27/8*(r+2)
    return (a*(3*r-2)+delta(r)**0.5/2)/(a*(r+2) + r**4), (a*(3*r-2)-delta(r)**0.5/2)/(a*(r+2) + r**4)

def proper_angle(r, cosa):
    num = 2*(1+cosa)/r + (1-3*cosa)
    den = (3-cosa) - 2*(1+cosa)/r
    return math.acos(-num/den)

radii = [ 0.1, 0.5, 1.0, 1.5, 2.0, 2.001, 2.999, 3.0, 4.0, 5.0, 6.0, 10.0 ]

for r in radii:
    cosa, cosa2 = cos_thresh(r)
    a01 =  math.acos(cosa)*180/math.pi
    a02 =  math.acos(cosa2)*180/math.pi
    try:
        limit_ang = proper_angle(r, cosa)*180/math.pi
        limit_ang2 = proper_angle(r, cosa2)*180/math.pi
    except:
        limit_ang = "NaN"
        limit_ang2 = "NaN"
    print("r =", r, "\tcos(a0) =", cosa, cosa2, "\ta0 =", a01, a02, "\ta0' =", limit_ang, limit_ang2)