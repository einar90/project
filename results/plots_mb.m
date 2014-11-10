clear all; close all; clc;

% Conversion factors
sm3PerDay2ccPerSec = 11.57;
mD2D = 0.001;
m2cm = 100.0;
bar2atm = 0.9869;

% Parameters
q = 4*150*sm3PerDay2ccPerSec; % cc/sec
mu = 0.5; % cP
k = 300.0*mD2D; % D
h = 30*m2cm; % cm

factor = (k*h)/(q*mu);

% Eclipse 10x10
close all
p_10x10 = load('eclipse/10x10-pressure.dat').*bar2atm;
p_0 = p_10x10(1,1); % Well-block pressure

p = (p_10x10-p_0).*factor;

x = [];
cc = 1 : 1 : 10;
rr1 = -1;
rr2 = 0;

for i=cc
  for j=cc
%     x(i,j) = sqrt((i-1)^2+(j-1)^2);
%     x(i,j) = sqrt((i+1)^2+(j+1)^2);
    x(i,j) = sqrt( ( i + rr1 )^2 + ( j + rr1 )^2 ) 
    + rr2;
  end
end

x

jj = 8;
pp = p(1 : jj, 1 : jj);
xx = x(1 : jj, 1 : jj);

pp = pp(:);
xx = ( xx(:) );
% Regression line for diagonal
% for i=1:4
%     x_diag(i) = sqrt((i-1)^2+(i-1)^2);
% end
% x_diag = xx(1 : ii);

ii = 20;
pfit = polyfit(log(xx(2:ii)), pp(2:ii), 1);
pval = polyval(pfit, log(xx(2:ii)));

hold on;
grid on
semilogx(xx(2:ii), pp(2:ii),'.')
% plot(log(xx(2:ii)), pp(2:ii),'.')

tta = .1:.1:6;
ttb = pfit(1).*log(tta) + pfit(2);
% semilogx(tta, ttb,'k')
plot(tta, ttb,'r')

-pfit(2)/pfit(1)
10^(-pfit(2)/pfit(1))

% axis([0 6 0 2.5]);

set(gca, 'XLim', [.1 10], 'YLim', [0 .6], ...
    'XScale', 'log', 'XTick', .1 : .1 : .4)



