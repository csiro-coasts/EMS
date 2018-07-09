[hed,a1]=hdrload('locs_z.ts');
[hed,b1]=hdrload('locs_s.ts');
[hed,a2]=hdrload('locb_z.ts');
[hed,b2]=hdrload('locb_s.ts');

close all;
subplot(3,1,1);
plot(a1(:,1),a1(:,2),'r-');
hold on;
plot(b1(:,1),b1(:,2),'b-');
ylabel('Elevation (m)');

subplot(3,1,2);
plot(a1(:,1),sqrt(a1(:,5).*a1(:,5)+a1(:,6).*a1(:,6)),'r-');
hold on;
plot(b1(:,1),sqrt(b1(:,5).*b1(:,5)+b1(:,6).*b1(:,6)),'b-');
ylabel('Surface speed (ms^-^1)');

subplot(3,1,3);
plot(a2(:,1),sqrt(a2(:,5).*a2(:,5)+a2(:,6).*a2(:,6)),'r-');
hold on;
plot(b2(:,1),sqrt(b2(:,5).*b2(:,5)+b2(:,6).*b2(:,6)),'b-');
ylabel('Bottom speed (ms^-^1)');
xlabel('Time (days)');
a=legend('z','sigma');

n = size(a1,1);
z_surf = sqrt(a1(n,5).*a1(n,5)+a1(n,6).*a1(n,6))
z_bot = sqrt(a2(n,5).*a2(n,5)+a2(n,6).*a2(n,6))
s_surf = sqrt(b1(n,5).*b1(n,5)+b1(n,6).*b1(n,6))
s_bot = sqrt(b2(n,5).*b2(n,5)+b2(n,6).*b2(n,6))

print -depsc -tiff test2.eps
