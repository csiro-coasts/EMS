
[hed,s1]=hdrload('loc1.ts');
[hed,s2]=hdrload('loc2.ts');
[hed,c1]=hdrload('loc1_quad.ts');
[hed,c2]=hdrload('loc2_quad.ts');
[hed,h1]=hdrload('loc1_hex.ts');
[hed,h2]=hdrload('loc2_hex.ts');
[hed,w1]=hdrload('loc1_quad2w.ts');
[hed,w2]=hdrload('loc2_quad2w.ts');
[hed,f1]=hdrload('loc1_hex2w.ts');
[hed,f2]=hdrload('loc2_hex2w.ts');


close all;
subplot(2,1,1);
plot(s1(:,1), s1(:,2), 'b-');
hold on;
plot(w1(:,1), w1(:,2), 'k-');
plot(c1(:,1), c1(:,2), 'r-');
plot(f1(:,1), f1(:,2), 'y-');
plot(h1(:,1), h1(:,2), 'g-');

ylabel('Eta (m)');

subplot(2,1,2);
plot(s1(:,1), s1(:,10), 'b-');
hold on;
plot(w1(:,1), w1(:,10), 'k-');
plot(c1(:,1), c1(:,10), 'r-');
plot(f1(:,1), f1(:,10), 'y-');
plot(h1(:,1), h1(:,10), 'g-');

xlabel('Time (days)');
ylabel('Temperature (^oC)');
title('Closed pool test');
a=legend('SHOC','Q2W','Quad','H2W','Hex','Location','SouthWest');

print -djpeg90 o1.jpg

