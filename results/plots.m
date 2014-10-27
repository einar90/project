clear all; close all; clc;

q = 1000; % bbl/day
mu = 0.96; % cP
k = 300.0e-3; % mD
h = 100.0; % ft
factor = q*mu / (k*h);

%% Eclipse 10x10
p_10x10 = load('eclipse/10x10-pressure.dat'); % psia
p_0 = p_10x10(1,1); % Well-block pressure

x = []; p = [];
for i=1:4
  for j=1:5
    x = [x sqrt(i^2+j^2)];
    p = [p (p_10x10(i,j)-p_0)/factor];
  end
end
    
semilogx(x,p,'.')
axis([0 6 0 2])
