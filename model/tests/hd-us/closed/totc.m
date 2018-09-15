[hed,aa]=hdrload('totals.ts');
close all;

subplot(5,1,1);
plot(aa(:,1),100*(aa(:,2) - aa(1,2))/aa(1,2));
title('Totals (% change)','FontSize',12);
ylabel('Mass (kg)');

subplot(5,1,2);
plot(aa(:,1),100*(aa(:,3) - aa(1,3))/aa(1,3));
ylabel('Volume (m^3)');

subplot(5,1,3);
plot(aa(:,1),100*(aa(:,4) - aa(1,4))/aa(1,4));
ylabel('Heat (^oCm^-^3)');

subplot(5,1,4);
plot(aa(:,1),100*(aa(:,5) - aa(1,5))/aa(1,5));
ylabel('Salt (psu.m^-^3)');
xlabel('Time (days)');

subplot(5,1,5);
plot(aa(:,1),100*(aa(:,6) - aa(1,6))/aa(1,6));
ylabel('Passive');
xlabel('Time (days)');

print -depsc o1.eps
print -djpeg90 o1.jpg
