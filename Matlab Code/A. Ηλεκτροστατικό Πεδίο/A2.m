clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                A.2 Parallel Plate Capacitor of Finite Width 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


platesWidth = 3 * 1e-2;           % width of plates in meters
platesThickness = 2 * 1e-3;       % thickness of plates in meters
platesDistance =  1 * 1e-2;       % distance between plates in meters
boxWidth =  5 * platesWidth;      % width of the computational box in meters
boxHeight = 5 * platesWidth;      % hight of the computational box in meters

refine = 1;                       % create denser mesh packed with triangles
jiggle = 'off';                   % jiggle mesh to spread nodes more evenly

[points, edges, triangles, geom] = ... 
geomplates(boxWidth, boxHeight, platesWidth, platesThickness, platesDistance, refine, jiggle);

numNodes = size(points,2);        % total number of nodes in mesh
numBoarder = size(edges,2);       % total number of edges at the boarder of mesh
numElements = size(triangles,2);  % total number of triangle elements



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


voltageUp = 50;                                  % voltage in the up plate
voltageDown = -50;                               % voltage in the down plate
InitValues = zeros(numNodes,1);                  % contains the Initial Potential Values at each node 
nodesKnown = zeros(numNodes,1);                  % contains 1 if the node has known value, else 0
nodesBoundary = union(edges(1, :), edges(2, :)); % take indexes of nodes that lie at the boundaries

isInsideBox =  (abs(points(1,nodesBoundary)) <= platesWidth / 2) & (abs(points(2,nodesBoundary)) <= platesThickness + platesDistance / 2);
isOutsideBox = (abs(points(1,nodesBoundary)) >= boxHeight / 2) | (abs(points(2,nodesBoundary)) >= boxWidth / 2);
isUp   = isInsideBox & (points(2,nodesBoundary) >=  platesDistance / 2);
isDown = isInsideBox & (points(2,nodesBoundary) <= -platesDistance / 2);

nodesUp = nodesBoundary(isUp);                   % take the nodes at the outter boundary
nodesDown = nodesBoundary(isDown);               % take the nodes at the inner boundary

InitValues(nodesUp) =   voltageUp;               % assign Dirichlet conditions at the outher boundary
InitValues(nodesDown) = voltageDown;             % assign Dirichlet conditions at the inner boundary
nodesKnown(union(nodesUp, nodesDown)) = 1;       % nodes at the boundaries have known values
nodesKnown(isOutsideBox) = 0;                    % assign zero valued Dirichlet (1) or Neuman (2) condition to outer boundary
numUnknowns = numNodes - sum(nodesKnown);         % degrees of freedom


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    PDE solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


e0 = 8.854 * 1e-12;               % vacuum permittivity
er = ones(1,numElements);         % relative permittivity outside plates

isbetweenPlates = (triangles(4, :) == 2);
er(isbetweenPlates) = 2.2;        % relative permittivity between plates
permittivity = er * e0;           % permittivity inside capacitor
method = 'direct';                % method of solver direct, bicg, gmres

[Potential, StiffnessGlobal] = pdesolver(points,triangles,permittivity,nodesKnown,InitValues,method);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PDE ploter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotMesh = 'off';
plotField = 'off'; 
plotPotential = 'off';

pdeploter(points, edges, triangles, geom, Potential, plotMesh, plotField, plotPotential);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


VoltDiff = voltageUp - voltageDown;                           % voltage difference in the 2 conductors
Energy = 1/2 * Potential' * (StiffnessGlobal * Potential);    % the approximated energy of the capacitor
CapacitancePredict = 2 * Energy / VoltDiff ^ 2;               % the approximated capacitance
CapacitanceTheory =  2.2 * e0 * platesWidth / platesDistance; % analytic calculation of capacitance

error = CapacitanceTheory - CapacitancePredict;          % error of aproximation
percent = abs(100 * error / CapacitanceTheory);           % percentege of error in relation to analytic capacitance
