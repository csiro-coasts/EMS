
[hed,s1]=hdrload('loc1.ts');
[hed,s2]=hdrload('loc2.ts');
[hed,c1]=hdrload('loc1_quad.ts');
[hed,c2]=hdrload('loc2_quad.ts');
[hed,h1]=hdrload('loc1_hex.ts');
[hed,h2]=hdrload('loc2_hex.ts');
[hed,w1]=hdrload('loc1_quad5w.ts');
[hed,w2]=hdrload('loc2_quad5w.ts');


close all;
subplot(2,1,1);
plot(s1(:,1), s1(:,2), 'b-');
hold on;
plot(c1(:,1), c1(:,2), 'k-');
plot(h1(:,1), h1(:,2), 'g-');
plot(w1(:,1), w1(:,2), 'r-');
ylabel('Eta (m)');

subplot(2,1,2);
plot(s1(:,1), s1(:,10), 'b-');
hold on;
plot(c1(:,1), c1(:,10), 'k-');
plot(h2(:,1), h2(:,10), 'g-');
plot(w1(:,1), w1(:,10), 'r-');
xlabel('Time (days)');
ylabel('Temperature (^oC)');
title('Closed basin test');
a=legend('SHOC','Quad','Hex','5W');

print -djpeg90 o1.jpg

