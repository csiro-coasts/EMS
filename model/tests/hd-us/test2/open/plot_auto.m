
% Plots an unstructured mesh using files <file>_e.txt and 
% <file>_c.txt output from shoc.

close all
load('in1_auto_e.txt');
x=in1_auto_e(:,1);
y=in1_auto_e(:,2);
line(x,y);
hold on

load('in1_auto_b.txt');
x=in1_auto_b(:,1);
y=in1_auto_b(:,2);
plot(x,y,'g-');

load('in1_auto_c.txt');
x=in1_auto_c(:,1);
y=in1_auto_c(:,2);
plot(x,y,'ro');


print -djpeg90 o1.jpg

