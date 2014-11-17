%% Define grid, rock and fluid data
nx = 20; ny = 20; nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G, 'Verbose', true);
rock.perm = repmat(100 .* milli*darcy, [G.cells.num, 1]);
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);
T = computeTrans(G, rock, 'Verbose', true);

%% Introduce wells
% The first well is vertical well (vertical is default):

% * completion in cells: cellsWell1
% * controlled by injection rate = 1.0  [m^3/d]
% * radius = 0.1.                       [m]
%
cellsWell1 =  1 : nx*ny : nx*ny*nz;
radius     = .1;
W = addWell([], G, rock, cellsWell1, 'Type', 'rate', ...
            'Val', 1.0/day(), 'Radius', radius, 'name', 'I');
disp('Well #1: '); display(W(1));

%%
% The second well is horizontal in the 'y' direction:
%
% * completion in cells: cellsWell2
% * controlled by bottom hole pressure, bhp = 1e5 [Pa]
% * radius = 0.1                                  [m]
%
cellsWell2 = nx : ny : nx*ny;
W = addWell(W, G, rock, cellsWell2, 'Type', 'bhp', ...
            'Val', 1.0e5, 'Radius', radius, 'Dir', 'y', 'name', 'P');
disp('Well #2: '); display(W(2));

%%
% We plot the wells to check if the wells are placed as we wanted them.
subplot(2,2,1), pos = get(gca,'Position'); clf
   plotGrid(G, 'FaceColor', 'none');
   view(3), camproj perspective, axis tight off, camlight headlight
   [ht, htxt, hs] = plotWell(G, W, 'radius', radius, 'height', 2);
   set(htxt, 'FontSize', 16);

%%
% Generate components of the linear system corresponding to the two wells, 
% initialize the solution structure (with correct bhp)
resSol = initState(G, W, 0);
display(resSol.wellSol);


%% Solve linear system
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off
resSol = incompTPFA(resSol, G, T, fluid, 'wells', W);

%% Report results
set(gca, 'Position', pos);  % move the current plot

subplot(2,2,2)
   plot(convertTo(-resSol.wellSol(2).flux, meter^3/day))
   title('Producer inflow profile [m^3/d]');

subplot(2,2,3)
   plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa));
   title('Pressure [bar]')
   view(3), camproj perspective, axis tight off, camlight headlight

subplot(2,2,4)
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   cf     = accumarray(cellNo, abs(resSol.flux(G.cells.faces(:,1))) );
   plotCellData(G, convertTo(cf, meter^3/day) );
   title('Sqrt of flux intensity [m^3/day]')
   view(3), camproj perspective, axis tight off, camlight headlight

%%
displayEndOfDemoMessage(mfilename)
