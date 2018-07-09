[hed,a1]=hdrload('loc_z1.ts');
[hed,b1]=hdrload('loc_z1.ts');
[hed,a2]=hdrload('loc_s2.ts');
[hed,b2]=hdrload('loc_s2.ts');

close all;
subplot(3,1,1);
plot(a1(:,1),a1(:,2),'r-');
hold on;
plot(b1(:,1),b1(:,2),'b-');
ylabel('West elevation (m)');

subplot(3,1,2);
plot(a2(:,1),a2(:,2),'r-');
hold on;
plot(b2(:,1),b2(:,2),'b-');
ylabel('East elevation (m)');

subplot(3,1,3);
plot(a1(:,1),sqrt(a2(:,3).*a2(:,3)+a2(:,4).*a2(:,4)),'r-');
hold on;
plot(b1(:,1),sqrt(b2(:,3).*b2(:,3)+b2(:,4).*b2(:,4)),'b-');
ylabel('Current speed (ms^-^1)');
xlabel('Time (days)');
a=legend('z','sigma');

n = size(a1,1);
dist = (9 * 2e3);
w_z = a1(n,2);
e_z = a2(n,2);
grad_z = (e_z - w_z) / dist
w_s = b1(n,2);
e_s = b2(n,2);
grad_s = (e_s - w_s) / dist

print -depsc -tiff test3.eps
