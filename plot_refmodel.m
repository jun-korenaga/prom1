% plot_refmodel.m
% 
% a sample script to show how to use calc_refmodel.m
%

tmin = 0;
tmax = 180;
dt = 1;
zmax = 300;
dz = 1;

[ts,zs,d,q,tt,zz,TT] = calc_refmodel(tmin,tmax,dt,zmax,dz);

figure(1);
subplot(3,1,1);
plot(ts,d);
axis ij;
xlabel('Age [Ma]');
ylabel('Depth [m]');

subplot(3,1,2);
plot(ts,q);
axis([tmin tmax 0 400]);
xlabel('Age [Ma]');
ylabel('Heat flow [mW/m^2]');

subplot(3,1,3); hold off;
pcolor(tt,zz,TT-273);
shading interp;
axis ij;
colorbar;
hold on;
[c,h]=contour(tt,zz,TT-273,[200:200:1000 1100:100:1400],'k-');
clabel(c,h);
xlabel('Age [Ma]');
ylabel('Depth [m]');
title('Temperature [^{\circ}C]');

% save to a file
tmp = [ts' d' q'];
save('ref_d_q.dat','tmp','-ascii');
tmp = [tt(:) zz(:) TT(:)];
save('ref_T.dat','tmp','-ascii');



 

