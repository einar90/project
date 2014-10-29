clear all; close all; clc;

% Conversion factors
sm3PerDay2ccPerSec = 11.57;
mD2D = 0.001;
m2cm = 100.0;
bar2atm = 0.9869;

% Parameters
q = 150*sm3PerDay2ccPerSec; % cc/sec
mu = 0.5; % cP
k = 300.0*mD2D; % D
h = 30*m2cm; % cm

factor = (k*h)/(q*mu);

%% Eclipse 10x10
p_10x10 = load('eclipse/10x10-pressure.dat').*bar2atm;
p_0 = p_10x10(1,1); % Well-block pressure

p = (p_10x10-p_0).*factor;

x = [];
for i=1:10
  for j=1:10
    x(i,j) = sqrt((i-1)^2+(j-1)^2);
  end
end

% Regression line for diagonal
for i=1:4
    x_diag(i) = sqrt((i-1)^2+(i-1)^2);
end
pfit = polyfit(x(2:4,2:4), p(2:4,2:4), 1);
pval = polyval(pfit, x_diag);
    

semilogx(x(2:4,2:4), p(2:4,2:4),'.')
hold on;
semilogx(x_diag,pval)
axis([0 6 0 2.5]);
