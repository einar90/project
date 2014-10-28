clear all; close all; clc;

% Conversion factors
bblDay2ccSec = 1.84;
ft2cm = 30.48;
psia2atm = 0.06805;

% Parameters
q = 1000*bblDay2ccSec; % cc/sec
mu = 0.96; % cP
k = 300.0e-3; % D
h = 100.0*ft2cm; % cm

factor = q*mu / (k*h);

%% Eclipse 10x10
p_10x10 = load('eclipse/10x10-pressure.dat')*psia2atm; % psia
p_0 = p_10x10(1,1); % Well-block pressure

p = (p_10x10-p_0)./factor;

x = [];
for i=1:10
  for j=1:10
    x(i,j) = sqrt(i^2+j^2);
    %p = [p (p_10x10(i,j)-p_0)/factor];
  end
end
    
%semilogx(x,p,'.')
plot(x(2,2), p(2,2),'.')
