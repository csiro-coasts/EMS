
[hed,s1]=hdrload('loc1.ts');
[hed,s2]=hdrload('loc2.ts');
[hed,c1]=hdrload('loc1_quad.ts');
[hed,c2]=hdrload('loc2_quad.ts');
[hed,h1]=hdrload('loc1_hex.ts');
[hed,h2]=hdrload('loc2_hex.ts');
[hed,w1]=hdrload('loc1_quad5w.ts');
[hed,w2]=hdrload('loc2_quad5w.ts');
[hed,f1]=hdrload('loc1_hex5w.ts');
[hed,f2]=hdrload('loc2_hex5w.ts');


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

ylabel('Salinity (PSU)');

%subplot(3,1,3);
%plot(s2(:,1), s1(:,10), 'b-');
%hold on;
%plot(c2(:,1), c1(:,10), 'k-');
%plot(h2(:,1), h2(:,10), 'g-');
%plot(w2(:,1), w2(:,10), 'r-');
%ylabel('Salinity (PSU)');

xlabel('Time (days)');
title('Estuarine test');
a=legend('SHOC','5WQ','Quad','5WH','Hex','Location','SouthEast');
%a=legend('SHOC','Quad','Hex','Location','SouthEast');

print -djpeg90 o1.jpg

