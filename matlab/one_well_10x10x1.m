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

%% Create 2D matrix from pressure solution vector and save it
P = reshape(convertTo(resSol.pressure, atm), dim, dim);
save([num2str(dim) 'x' num2str(dim) '-pressure.dat'], 'P', '-ascii')
