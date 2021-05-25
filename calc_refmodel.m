% calc_refmodel.m
% 
% A Matlab code for calculating the reference model 
% of Korenaga et al. (JGR, 2021)
%
% input parameters:
% - t_start, t_end, dt: time range in million years and time increment
% - zmax, dz: maximum depth and depth increment
%
% output parameters:
% - ts: time
% - zs: depth
% - d: seafloor depth [m]
% - q: surface heat flow [mW/m^2]
% - tt: time grid
% - zz: depth grid
% - TT: temperature [K]

function [ts,zs,d,q,tt,zz,TT] = calc_refmodel(t_start,t_end,dt,zmax,dz)

ts = t_start:dt:t_end;
if (ts(1) == 0)
  ts(1) = 1e-6; % to avoid singularity
end

sqt = sqrt(ts);
t32 = ts.^(3/2);
t2 = ts.^2;

zs = 0:dz:zmax;
[tt,zz] = meshgrid(ts,zs);

%
% seafloor depth
%
d0 = 2600;
d1 = 409; f_TC = 0.85; d2 = 930; e2 = 0.018;
p1 = 32.85; p2 = -18.39; p3 = 0.3023; p4 = -0.0054;

d = d0 + d1*f_TC*sqt + d2*tanh(e2*ts) ...
    + p1*sqt + p2*ts + p3*t32 + p4*t2;

%
% surface heat flow
%
c0 = 338.4; c1 = 66.7; c2 = -8.26; c3 = 0.53; c4 = -0.013;


H = 2.3e-12;
gamma = 3.154e13;
Cp = 1200;
DT = 1350;
dT = 0.1*ts; % secular cooling of 100 K/Gyr

q = 1./sqt.*(1+2*H*gamma*ts/(Cp*DT)+dT/DT).*(c0+c1*sqt+c2*ts+c3*t32+c4*t2);

%
% thermal structure
%
Ts = 273;
a1 = 0.602e-3; a2 = -6.045e-10;
kappa7 = 3.45e-7;
kappa0 = 2.23e-6; zref = 1e5;
b = [-1.255 9.944 -25.0619 32.2944 -22.2017 7.7336 -1.0622];

zzm = zz*1e3; % now in m
tts = tt*gamma; % now in s

sum_b = 0;
for i=0:6
  sum_b = sum_b + b(i+1)*(zzm/zref).^(i/2);
end
kappaz = kappa0*sum_b;
kappaz(zzm<7e3) = kappa7;

T_KK16 = Ts + DT*erf(zzm./(2*sqrt(kappaz.*tts))) ...
	 + a1*zzm + a2*zzm.^2;

kappa = 1e-6;
fac1 = zzm./(2*sqrt(kappa*tts));
fac2 = zzm.^2/(2*kappa);
erf_fac1 = erf(fac1);

T1 = Ts + DT*erf_fac1;
T2 = T1 + H/Cp*((tts+fac2).*erf_fac1 ...
		+ zzm.*sqrt(tts/(pi*kappa)).*exp(-fac1.^2) ...
		- fac2) ...
     + dT.*erf_fac1;

TT = T_KK16.*T2./T1;
