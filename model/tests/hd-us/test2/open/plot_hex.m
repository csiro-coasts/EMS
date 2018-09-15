
% Plots an unstructured mesh using files <file>_e.txt and 
% <file>_c.txt output from shoc.

close all
load('in1_hex_e.txt');
x=in1_hex_e(:,1);
y=in1_hex_e(:,2);
line(x,y);
hold on;

load('in1_hex_c.txt');
x=in1_hex_c(:,1);
y=in1_hex_c(:,2);
plot(x,y,'ro');

print -djpeg90 o1.jpg

