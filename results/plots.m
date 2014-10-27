clear all; close all; clc;

q = 1000; % bbl/day
mu = 0.96; % cp
k = 300.0e-3; % F
h = 100.0; % ft

%% Eclipse 10x10
p_10x10 = load('eclipse/10x10-pressure.dat');
p_0 = p_10x10(1,1); % Well-block pressure

% Diagonal pressures for plotting
for i=2:3
    p_diag(i) = (p_10x10(i,i)-p_0)/(q*mu/k/h); 
    p_xedge(i) = (p_10x10(1,i)-p_0)/(q*mu/k/h); 
    p_yedge(i) = (p_10x10(i,1)-p_0)/(q*mu/k/h);
    x(i) = sqrt(2*i^2);
end

semilogx(x,p_diag,'.')
hold on;
semilogx(x,p_xedge,'x')
axis([0 6 0 .6])
