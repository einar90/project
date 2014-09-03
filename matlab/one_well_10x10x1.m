clear all; clc;
%% Define grid, rock and fluid data
dim = 20;
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
P = reshape(resSol.pressure, dim, dim);

%% Pressures on diagonal from producer towards center
for i = 2:dim/2
  diagP(i-1) = P(i,i);
  r(i-1) = sqrt(2*(i^2));
end

%% Regression line
pfit = polyfit(r, diagP, 1);
pval = polyval(pfit, r);

%% Report results
subplot(1,2,1)
  surf(P)
  title('Pressure [bar]')
  view(2), axis tight, colorbar
  
subplot(1,2,2)
  semilogx(r, diagP, '*')
  title('Pressure vs. radius')
  xlabel('$$r/\Delta x = \sqrt{i^2 + j^2}$$','interpreter','latex')
  ylabel('$$P$$' ,'interpreter','latex')
  ylim([0 max(diagP)]);
  xlim([1, max(r)])
