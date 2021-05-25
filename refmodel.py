# calc_refmodel.m
# 
# A Phython module for calculating the reference model 
# of Korenaga et al. (JGR, 2021)
#
# input parameters:
# - t_start, t_end, dt: time range in million years and time increment
# - zmax, dz: maximum depth and depth increment
#
# output parameters:
# - ts: time
# - zs: depth
# - d: seafloor depth [m]
# - q: surface heat flow [mW/m^2]
# - tt: time grid
# - zz: depth grid
# - TT: temperature [K]

def calc(t_start,t_end,dt,zmax,dz):
    """
    """
    import numpy as np
    import scipy.special as ss
    import math

    ts = np.arange(t_start,t_end+dt,dt+0.0)
    if ts[0] == 0:
        ts[0] = 1e-6 # to avoid singularity

    sqt = np.sqrt(ts)
    t32 = ts**(3/2)
    t2 = ts**2

    zs = np.arange(0,zmax+dz,dz+0.0)
    tt, zz = np.meshgrid(ts,zs)

    #
    # seafloor depth
    #
    d0 = 2600;
    d1, f_TC, d2, e2 =  409, 0.85, 930, 0.018
    p1, p2, p3, p4 = 32.85, -18.39, 0.3023, -0.0054

    d = d0 + d1*f_TC*sqt + d2*np.tanh(e2*ts) \
        + p1*sqt + p2*ts + p3*t32 + p4*t2

    #
    # surface heat flow
    #
    c0, c1, c2, c3, c4 = 338.4, 66.7, -8.26, 0.53, -0.013
    H, gamma, Cp, DT = 2.3e-12, 3.154e13, 1200, 1350
    dT = 0.1*ts # secular cooling of 100 K/Gyr

    q = (1/sqt)*(1+2*H*gamma*ts/(Cp*DT)+dT/DT)*(c0+c1*sqt+c2*ts+c3*t32+c4*t2)

    #
    # thermal structure
    #
    Ts = 273
    a1, a2 = 0.602e-3, -6.045e-10
    kappa7 = 3.45e-7
    kappa0, zref = 2.23e-6, 1e5
    b = np.array([-1.255, 9.944, -25.0619, 32.2944, -22.2017, 7.7336, -1.0622])

    zzm = zz*1e3 # now in m
    tts = tt*gamma # now in s

    sum_b = np.zeros(zzm.shape);
    for i in np.arange(0,7,1):
        sum_b = sum_b + b[i]*(zzm/zref)**(i/2);
    kappaz = kappa0*sum_b;
    kappaz[zzm<7e3] = kappa7;

    T_KK16 = Ts + DT*ss.erf(zzm/(2*np.sqrt(kappaz*tts))) \
	     + a1*zzm + a2*zzm**2;

    kappa = 1e-6
    fac1 = zzm/(2*np.sqrt(kappa*tts))
    fac2 = zzm**2/(2*kappa)
    erf_fac1 = ss.erf(fac1);

    T1 = Ts + DT*erf_fac1;
    T2 = T1 + H/Cp*((tts+fac2)*erf_fac1 \
		+ zzm*np.sqrt(tts/(math.pi*kappa))*np.exp(-fac1**2) \
		- fac2) \
                + dT*erf_fac1;

    TT = T_KK16*T2/T1;

    return ts, d, q, tt, zz, TT


