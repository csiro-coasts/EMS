
% Plots an unstructured mesh using files <file>_e.txt and 
% <file>_c.txt output from shoc.

close all
load('est_hex_e.txt');
plot(est_hex_e(:,1),est_hex_e(:,2),'r-');
hold on;

load('est_hex_c.txt');
plot(est_hex_c(:,1),est_hex_c(:,2),'b.');

load('boundary.txt');
plot(boundary(:,1),boundary(:,2),'g-');

print -djpeg90 o1.jpg

