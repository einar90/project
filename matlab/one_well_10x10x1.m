clear all; clc;
%% Define grid, rock and fluid data
dim = 10;
G = cartGrid([dim, dim, 1]);
G = computeGeometry(G, 'Verbose', true);
rock.perm = repmat(100 .* milli*darcy, [G.cells.num, 1]);
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);
T = computeTrans(G, rock, 'Verbose', true);

%% Wells
r = .1;
W = [];

%% Production well
cellsWell1 =  [1];
W = addWell(W, G, rock, cellsWell1, 'Type', 'rate', ...
            'Val', -1.0/day(), 'Radius', r, 'name', 'P');

%% Injection well
cellsWell2 = [dim^2];
W = addWell(W, G, rock, cellsWell2, 'Type', 'rate', ...
            'Val', 1.0/day(), 'Radius', r, 'name', 'I');


%% Generate linear system, initialize solution structure
resSol = initState(G, W, 0);


%% Solve linear system
gravity off
resSol = incompTPFA(resSol, G, T, fluid, 'wells', W);

%% Create 2D matrix from pressure solution vector
P = reshape(convertTo(resSol.pressure, atm), dim, dim);

%% Pressures on diagonal from producer towards center
for i = 2:dim/2
  diagP(i-1) = P(i,i);
  r(i-1) = sqrt(2*(i^2));
end

%% Pressure line for plotting
sm3PerDay2ccPerSec = 11.57;
k = 100*0.001;
h = 100;
q = 1*sm3PerDay2ccPerSec;
mu = 1;
factor = k*h/(q*mu);
p_plot = (P - P(1,1)).*factor;

for i=1:dim
    for j=1:dim
        x(i,j) = sqrt((i-1)^2 + (j-1)^2);
    end
end

% Cut matrices
p_plot = p_plot(2:3,2:3);
x = x(2:3,2:3);

%% Regression line
pfit = polyfit(x, p_plot, 1);
pval = polyval(pfit, x);

%% Report results
subplot(1,2,1)
  surf(P)
  title('Pressure [bar]')
  view(2), axis tight, colorbar
  
subplot(1,2,2)
  plot(x, p_plot, '*')
  title('Pressure vs. radius')
  xlabel('$$r/\Delta x = \sqrt{i^2 + j^2}$$','interpreter','latex')
  ylabel('$$P$$' ,'interpreter','latex')
  ylim([0 max(diagP)]);
  xlim([0, 6])
