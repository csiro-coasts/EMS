[hed,aa]=hdrload('totals.ts');
close all;

subplot(2,1,2);
plot(aa(:,1),100*(aa(:,6) - aa(1,6))/aa(1,6));
ylabel('Passive');
xlabel('Time (days)');

print -depsc o1.eps
print -djpeg90 o1.jpg
