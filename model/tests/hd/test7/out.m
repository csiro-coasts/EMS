[hed,aa]=hdrload('loc_z.ts');
[hed,bb]=hdrload('loc_s.ts');
[hed,cc]=hdrload('loc_e.ts');

close all;
subplot(3,1,1);
plot(aa(:,1),aa(:,2),'r-');
hold on;
plot(bb(:,1),bb(:,2),'b-');
plot(cc(:,1),cc(:,2),'g-');
ylabel('Elevation (m)');

subplot(3,1,2);
plot(aa(:,1),aa(:,3),'r-');
hold on;
plot(bb(:,1),bb(:,3),'b-');
plot(cc(:,1),cc(:,3),'g-');
ylabel('U1 (ms^-^1)');

subplot(3,1,3);
plot(aa(:,1),aa(:,4),'r-');
hold on;
plot(bb(:,1),bb(:,4),'b-');
plot(cc(:,1),cc(:,4),'g-');
ylabel('U2 (ms^-^1)');
xlabel('Time (days)');
a=legend('z','sigma','Mapped','Location','NorthWest');

print -depsc -tiff test7.eps

%%%%%%%%%%%%%%%%%%%%%%%%%
% Totals
kpro=input('Enter any key to continue : '); 
[hed,aa]=hdrload('totals.ts');
close all;

subplot(4,1,1);
plot(aa(:,1),100*(aa(:,2) - aa(1,2))/aa(1,2));
title('Totals (% change)','FontSize',12); 
ylabel('Mass (kg)');

subplot(4,1,2);
plot(aa(:,1),100*(aa(:,3) - aa(1,3))/aa(1,3));
ylabel('Volume (m^3)');

subplot(4,1,3);
plot(aa(:,1),100*(aa(:,4) - aa(1,4))/aa(1,4));
ylabel('Heat (deg C m^3)');

subplot(4,1,4);
plot(aa(:,1),100*(aa(:,5) - aa(1,5))/aa(1,5));
ylabel('Salt (psu m^3)');
xlabel('Time (days)');

print -depsc -tiff total.eps

kpro=input('Enter any key to continue : '); 

subplot(3,1,1);
plot(aa(:,1),100*(aa(:,6) - aa(1,6))/aa(1,6));
title('Totals (% change)','FontSize',12);
ylabel('Sand (kg)');

subplot(3,1,2);
plot(aa(:,1),100*(aa(:,7) - aa(1,7))/aa(1,7));
ylabel('Silt (kg)');

subplot(3,1,3);
plot(aa(:,1),100*(aa(:,8) - aa(1,8))/aa(1,8));
ylabel('Fine (kg)');
xlabel('Time (days)');

print -depsc -tiff sediment.eps

%%%%%%%%%%%%%%%%%%%%%%%%%
% Flushing
kpro=input('Enter any key to continue : ');
[hed,aa]=hdrload('flushing.ts');
close all;

n=3;       % Order of polynomial fit
t=13.06;   % Actual flushing time
subplot(2,1,1);
[p,s]=polyfit(aa(:,1),aa(:,3),n);
n=aa(size(aa(:,1),1),1);
m=n/100.0;
b=[0:m:n];
bb=polyval(p,b);
hold off;
plot(aa(:,1),aa(:,3));
hold on;
plot(b,bb,'r-');
e=2.718281828;
y=[1/e 1/e];
x=[0 n];
plot(aa(:,1),exp(-aa(:,1)/t),'g-');
%plot(x,y,'k:');
%axis([0,1,0,1]);
xlabel('Time (days)','FontSize',12);
ylabel('Normalized Total Mass','FontSize',12);
%a=legend('Data','Polynomial fit','Exponential fit','e-folding','SouthEastOutside');
a=legend('Data','Polynomial fit','Exponential fit',0);
print -depsc -tiff flush.eps

%%%%%%%%%%%%%%%%%%%%%%%%%
% Sourcesink
kpro=input('Enter any key to continue : '); 
[hed,aa]=hdrload('totals_diag.ts');
close all;

r=100;
t=[0:.1:10];
m=r*86400*t;

subplot(2,1,1);
plot(aa(:,1),aa(:,6),'r-');
hold on;
plot(t,m,'k--');
title('Mass balance of the sourcesink','FontSize',12);
xlabel('Time (days)','FontSize',12);
ylabel('Total Mass','FontSize',12);
a=legend('Model','Theory','SouthEast');

m1=r*86400*10;
sz=size(aa,1);
v1=aa(sz,6);
percent_mass_error=100-100*m1/v1

%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple windows
kpro=input('Enter any key to continue : ');
[hed,aa]=hdrload('loc_1w.ts');
[hed,bb]=hdrload('loc_4w.ts');
close all;

subplot(4,1,1);
plot(aa(:,1),aa(:,2),'r-');
hold on;
plot(bb(:,1),bb(:,2),'b--');
ylabel('Eta (m)');

subplot(4,1,2);
plot(aa(:,1),aa(:,3),'r-');
hold on;
plot(aa(:,1),aa(:,4),'k-');
plot(bb(:,1),bb(:,3),'b--');
plot(bb(:,1),bb(:,4),'g--');
ylabel('2D Vel (ms^-^1)');

subplot(4,1,3);
plot(aa(:,1),aa(:,5),'r-');
hold on;
plot(aa(:,1),aa(:,6),'k-');
plot(bb(:,1),bb(:,5),'b--');
plot(bb(:,1),bb(:,6),'g--');
ylabel('3D Vel (ms^-^1)');

subplot(4,1,4);
plot(aa(:,1),aa(:,11),'r-');
hold on;
plot(bb(:,1),bb(:,11),'b--');
ylabel('Temp (^oC)');
xlabel('Time (days)');

print -depsc -tiff windows.eps
