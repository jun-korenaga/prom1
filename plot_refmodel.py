# plot_refmodel.py
#
# a sample script to show how to use refmodel.py

import numpy as np, matplotlib.pyplot as plt
import refmodel as refm

tmin, tmax, dt = 0, 180, 1
zmax, dz = 300, 1
ts, d, q, tt, zz, TT = refm.calc(tmin,tmax,dt,zmax,dz)

plt.figure();
plt.plot(ts,d)
plt.xlabel('Age [Ma]'); plt.ylabel('Depth [m]');
plt.axis([tmin, tmax, 2500, 6000])
plt.gca().invert_yaxis()

plt.figure();
plt.plot(ts,q)
plt.xlabel('Age [Ma]'); plt.ylabel('Heat flow [mW/m$^2$]');
plt.axis([tmin, tmax, 0, 400])

plt.figure();
plt.pcolormesh(tt, zz, TT, shading='gouraud')
CS = plt.contour(tt, zz, TT, \
                 [200, 400, 600, 800, 1000, 1100, 1200, 1300, 1400])
plt.clabel(CS, inline=True, fontsize=10, fmt='%d')
plt.xlabel('Age [Ma]'); plt.ylabel('Depth [m]');
plt.gca().invert_yaxis()

# save to a file
outf = open('pref_d_q.dat','w')
for i in np.arange(0,ts.size,1):
    outf.write("{:f}\t{:f}\t{:f}\n".format(ts[i],d[i],q[i]))
outf.close()

nt, nz = tt.shape;

outf = open('pref_T.dat','w')
for j in np.arange(0,nz,1):
    for i in np.arange(0,nt,1):
        outf.write("{:f}\t{:f}\t{:f}\n".format(tt[i,j],zz[i,j],TT[i,j]))
outf.close()
    

