clear all; clc; close all;
%% Define grid, rock and fluid data
dim = 53;
x = zeros(1, dim+1);
x(1:3) = [0 3.75 15];
for i=4:dim-1
  x(i) = x(i-1)+30;
end
x(end-1) = x(end-2)+11.25;
x(end)   = x(end-1)+3.75;

G = tensorGrid(x, x, [0 30]); % dim,dim,1 blocks, 30x30x30m each
plotGrid(G); axis([0 300 0 300]);
G = computeGeometry(G, 'Verbose', true);
rock.perm = repmat(1 .*darcy, [G.cells.num, 1]);
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

%% Create 2D matrix from pressure solution vector and save it
P = reshape(convertTo(resSol.pressure, atm), dim, dim);
save([num2str(dim) 'x' num2str(dim) '-pressure.dat'], 'P', '-ascii')

