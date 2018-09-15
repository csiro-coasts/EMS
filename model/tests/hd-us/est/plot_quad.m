
% Plots an unstructured mesh using files <file>_e.txt and 
% <file>_c.txt output from shoc.

close all
load('est_quad_e.txt');
plot(est_quad_e(:,1),est_quad_e(:,2),'b-');
hold on;

load('est_quad_c.txt');
plot(est_quad_c(:,1),est_quad_c(:,2),'ro');

%load('est_quad_b.txt');
%plot(est_quad_b(:,1),est_quad_b(:,2),'g-');

print -djpeg90 o1.jpg

