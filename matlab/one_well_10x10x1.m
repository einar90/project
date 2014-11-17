clear all; clc; close all;
%% Define grid, rock and fluid data
dim = 11;
G = cartGrid([dim, dim, 1], [30, 30, 30]); % dim,dim,1 blocks, 30x30x30m each
G = computeGeometry(G, 'Verbose', true);
rock.perm = repmat(.3 .*darcy, [G.cells.num, 1]);
fluid     = initSingleFluid('mu' ,    0.5*centi*poise     , ...
                            'rho', 1000*kilogram/meter^3);
T = computeTrans(G, rock, 'Verbose', true);

%% Wells
r = .3;
W = [];

%% Production well
cellsWell1 =  [1];
W = addWell(W, G, rock, cellsWell1, 'Type', 'rate', ...
            'Val', -150.0/day(), 'Radius', r, 'name', 'P');

%% Injection well
cellsWell2 = [dim^2];
W = addWell(W, G, rock, cellsWell2, 'Type', 'rate', ...
            'Val', 150.0/day(), 'Radius', r, 'name', 'I');


%% Generate linear system, initialize solution structure
resSol = initState(G, W, 0);


%% Solve linear system
gravity off
resSol = incompTPFA(resSol, G, T, fluid, 'wells', W);

%% Create 2D matrix from pressure solution vector
P = reshape(convertTo(resSol.pressure, atm), dim, dim);
save([num2str(dim) 'x' num2str(dim) '-pressure.dat'], 'P', '-ascii')


%% Pressure line for plotting
sm3PerDay2ccPerSec = 11.57;
k = .3; % D
h = 30.0*100; % cm
q = 150.0*sm3PerDay2ccPerSec; % cc/sec
mu = .5; % cP
factor = k*h/(q*mu);
p_plot = (P - P(1,1)).*factor;

for i=1:dim
    for j=1:dim
        x(i,j) = sqrt((i-1)^2 + (j-1)^2);
    end
end



%% Regression line
xx = x(1:floor(dim/2),1:floor(dim/2));
pp = p_plot(floor(1:dim/2),floor(1:dim/2));
xx = xx(:);
pp = pp(:);
pfit = polyfit(log(xx(2:end)), pp(2:end), 1);
pval = polyval(pfit, xx);
x_reg = [.1:.01:10];

%% Report results
subplot(1,2,1)
  surf(P)
  title('Pressure [bar]')
  view(2), axis tight, colorbar
  
subplot(1,2,2)
  semilogx(xx, pp, '*')
  title('Pressure vs. radius')
  xlabel('$$r/\Delta x = \sqrt{i^2 + j^2}$$','interpreter','latex')
  ylabel('$$P$$' ,'interpreter','latex')
  ylim([0 5]);
  xlim([10e-2, .6e1])
  hold on 
  semilogx(x_reg, pfit(2)+pfit(1).*log(x_reg))
