
% Plots an unstructured mesh using files <file>_e.txt and 
% <file>_c.txt output from compas.

close all

% Plot edges
load('auto_e.txt');
plot(auto_e(:,1),auto_e(:,2),'b-');
hold on;

% Plot cell centres
load('auto_c.txt');
plot(auto_c(:,1),auto_c(:,2),'r.');

% Include original coastline
load('setas.cst');
%plot(setas(:,1),setas(:,2),'k-');

% Include processed coastline
load('cetas_p_out.txt');
plot(cetas_p_out(:,1),cetas_p_out(:,2),'k-');

