[hed,aa]=hdrload('loc_z.ts');
[hed,bb]=hdrload('loc_msl_z.ts');

close all;
subplot(3,1,1);
plot(aa(:,1),aa(:,2),'r-');
hold on;
plot(bb(:,1),bb(:,2),'b-');
ylabel('Elevation (m)');

subplot(3,1,2);
plot(aa(:,1),aa(:,3),'r-');
hold on;
plot(bb(:,1),bb(:,3),'b-');
ylabel('U1 (ms^-^1)');

subplot(3,1,3);
plot(aa(:,1),aa(:,4),'r-');
hold on;
plot(bb(:,1),bb(:,4),'b-');
ylabel('U2 (ms^-^1)');
xlabel('Time (days)');
a=legend('below msl','above msl');

print -depsc -tiff test8.eps
