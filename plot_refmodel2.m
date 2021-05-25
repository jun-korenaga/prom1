% plot_refmodel2.m
% 
% another sample script to show how to use calc_refmodel.m
% 
% this version compares with observations and creates
% something like Figures 7a and 8 of Korenaga et al. (2021).
%

tmin = 0;
tmax = 180;
dt = 1;
zmax = 300;
dz = 1;

[ts,zs,d,q,tt,zz,TT] = calc_refmodel(tmin,tmax,dt,zmax,dz);

% load age-depth data of normal seafloor
data = load('normal_sa_depth_hist.dat');
sgrid = data(:,1);
dgrid = data(:,2);
pgrid = data(:,3);
ns = length(unique(sort(sgrid)));
nd = length(unique(sort(dgrid)));
sgrid = reshape(sgrid,nd,ns);
dgrid = reshape(dgrid,nd,ns);
pgrid = reshape(pgrid,nd,ns);

% load age-heatflow data of normal seafloor
data = load('hf_quartile_filHFnormal.dat');
range = logical(~isnan(data(:,3)) & data(:,2)>4);
t0 = data(range,1);
t1 = t0+2.5;
q2 = data(range,3);
q1 = data(range,4);
q3 = data(range,5);
LNdata = load('lister_nagihara.dat');

figure(2);
subplot(2,1,1); hold off;
pcolor(sgrid,dgrid,pgrid); axis ij; 
shading interp; colorbar;
axis([0 13.5 2000 7000]); hold on;

plot(sqrt(ts),d,'w-','LineWidth',2);
xlabel('Age^{1/2} [Ma^{1/2}]');
ylabel('Depth [m]');

subplot(2,1,2); hold off;
plot(ts,q,'k','LineWidth',2);
axis([tmin tmax 0 150]);
xlabel('Age [Ma]');
ylabel('Heat flow [mW/m^2]');
hold on;
errorbar(LNdata(:,1),LNdata(:,2),LNdata(:,3),'r^');
for i=1:length(t0)
  plot([t0(i) t0(i) t1(i) t1(i) t0(i)],...
       [q1(i) q3(i) q3(i) q1(i) q1(i)],'b-');
  plot(t0(i)+1.25,q2(i),'co');
end
legend('KK21', 'Lister+Nagihara','IQR','median');




 

