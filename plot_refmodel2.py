# plot_refmodel2.py
#
# a sample script to show how to use refmodel.py
# with some data to compare

import numpy as np, matplotlib.pyplot as plt
import refmodel as refm

tmin, tmax, dt = 0, 180, 1
zmax, dz = 300, 1
ts, d, q, tt, zz, TT = refm.calc(tmin,tmax,dt,zmax,dz)

# load age-depth data of normal seafloor
data = np.loadtxt('normal_sa_depth_hist.dat')
sgrid = data[:,0]
dgrid = data[:,1]
pgrid = data[:,2]
ns = np.unique(np.sort(sgrid)).size
nd = np.unique(np.sort(dgrid)).size
sgrid = sgrid.reshape((ns,nd))
dgrid = dgrid.reshape((ns,nd))
pgrid = pgrid.reshape((ns,nd))

# load age-heatflow data of normal seafloor
data = np.loadtxt('hf_quartile_filHFnormal.dat');
range = ((~np.isnan(data[:,2])) & (data[:,1]>4))
t0 = data[range,0]
t1 = t0+2.5
q2 = data[range,2]
q1 = data[range,3]
q3 = data[range,4]
LNdata = np.loadtxt('lister_nagihara.dat')

plt.figure();
plt.pcolormesh(sgrid, dgrid, pgrid, shading='gouraud')
plt.plot(np.sqrt(ts),d,'w-')
plt.xlabel('Age$^{1/2}$ [Ma$^{1/2}$]'); plt.ylabel('Depth [m]');
plt.axis([0, 13.5, 2000, 7000])
plt.gca().invert_yaxis()

plt.figure();
plt.plot(ts,q,'k-', label='PROM1')
plt.xlabel('Age [Ma]'); plt.ylabel('Heat flow [mW/m$^2$]');
plt.axis([tmin, tmax, 0, 150])
for i in np.arange(0,t0.size,1):
    line, = plt.plot([t0[i], t0[i], t1[i], t1[i], t0[i]], \
                     [q1[i], q3[i], q3[i], q1[i], q1[i]], \
                     'b-', linewidth=0.5)
    if i==1:
      line.set_label('IQR')
    line2, = plt.plot(t0[i]+1.25,q2[i],'co',markersize=2)
    if i==1:
      line2.set_label('median')

plt.errorbar(LNdata[:,0],LNdata[:,1],LNdata[:,2],\
             np.zeros(LNdata[:,2].shape), 'none', \
             ecolor='red', label='Lister+Nagihara')
plt.legend()
