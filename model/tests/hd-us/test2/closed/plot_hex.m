
% Plots an unstructured mesh using files <file>_e.txt and 
% <file>_c.txt output from shoc.
load('in1_hex_e.txt');
load('in1_hex_c.txt');
a1=in1_hex_e;
b1=in1_hex_c;

close all

plot(a1(:,1),a1(:,2),'b-');
hold on;

plot(b1(:,1),b1(:,2),'ro');

print -djpeg90 o1.jpg

