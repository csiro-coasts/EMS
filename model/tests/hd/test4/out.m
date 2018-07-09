[hed,aa]=hdrload('loc_z.ts');
[hed,bb]=hdrload('loc_s.ts');

close all;
subplot(3,1,1);
plot(aa(:,1),aa(:,2),'r-');
hold on;
plot(bb(:,1),bb(:,2),'b-');
ylabel('Elevation (m)');

d=52.5;
r=0.0005;
tsx=0.1;
rho=1024.76;
u=(d*tsx*(1-exp(-aa(:,1)*86400*r/d))/(rho*r))/d;

sz = size(aa,1);
eta_z = aa(sz,2)
U_z = aa(sz,3)
V_z = aa(sz,4)
subplot(3,1,2);
plot(aa(:,1),aa(:,3),'r-');
hold on;
plot(bb(:,1),bb(:,3),'b-');
plot(aa(:,1),u,'g-');
ylabel('U1 (ms^-^1)');

sz = size(bb,1);
eta_s = bb(sz,2)
U_s = bb(sz,3)
V_s = bb(sz,4)
subplot(3,1,3);
plot(aa(:,1),aa(:,4),'r-');
hold on;
plot(bb(:,1),bb(:,4),'b-');
ylabel('U2 (ms^-^1)');
xlabel('Time (days)');
a=legend('z','sigma');

print -depsc -tiff test4.eps
