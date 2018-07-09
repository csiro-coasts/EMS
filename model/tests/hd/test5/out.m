[hed,a1]=hdrload('loc1_z.ts');
[hed,a2]=hdrload('loc2_z.ts');
[hed,a3]=hdrload('loc3_z.ts');
[hed,a4]=hdrload('loc4_z.ts');
[hed,a5]=hdrload('loc5_z.ts');
[hed,a6]=hdrload('loc6_z.ts');
[hed,a7]=hdrload('loc7_z.ts');
[hed,a8]=hdrload('loc8_z.ts');
[hed,a9]=hdrload('loc9_z.ts');
[hed,a10]=hdrload('loc10_z.ts');

[hed,b1]=hdrload('loc1_s.ts');
[hed,b2]=hdrload('loc2_s.ts');
[hed,b3]=hdrload('loc3_s.ts');
[hed,b4]=hdrload('loc4_s.ts');
[hed,b5]=hdrload('loc5_s.ts');
[hed,b6]=hdrload('loc6_s.ts');
[hed,b7]=hdrload('loc7_s.ts');
[hed,b8]=hdrload('loc8_s.ts');
[hed,b9]=hdrload('loc9_s.ts');
[hed,b10]=hdrload('loc10_s.ts');

close all;
n = size(a1,1);
dist = [5 15 25 35 45 55 65 75 85 95];
d1 = [.0301 0.0189 0.0137 0.0103 0.0077 0.0057 0.0039 0.0025 0.0012 0.0001];
d2 = [a1(n,2) a2(n,2) a3(n,2) a4(n,2) a5(n,2) a6(n,2) a7(n,2) a8(n,2) a9(n,2) a10(n,2)];
d3 = [b1(n,2) b2(n,2) b3(n,2) b4(n,2) b5(n,2) b6(n,2) b7(n,2) b8(n,2) b9(n,2) b10(n,2)];

subplot(3,1,1);
plot(dist,d1,'r-');
hold on;
plot(dist,d2,'b-');
plot(dist,d3,'g-');
a=legend('theory','z','sigma');
ylabel('Elevation (m)');
xlabel('Distance (km)');

sz = size(a5,1);
U_z = a5(sz,3)
U_s = b5(sz,3)
subplot(3,1,2);
plot(a5(:,1),a5(:,3),'r-');
hold on;
plot(b5(:,1),b5(:,3),'b-');
ylabel('U1 (ms^-^1)');

sz = size(a5,1);
V_z = a5(sz,4)
V_s = b5(sz,4)
subplot(3,1,3);
plot(a5(:,1),a5(:,4),'r-');
hold on;
plot(b5(:,1),b5(:,4),'b-');
ylabel('U2 (ms^-^1)');
xlabel('Time (days)');
a=legend('z','sigma');

print -depsc -tiff test5.eps
