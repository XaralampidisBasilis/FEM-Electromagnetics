clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   A.1 Coaxial Cable with Air fielectric 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Z0 = 50;                                 % characteristic resistance of the coaxial cable  
radiusOut = 3.5 / 2 * 1e-3;              % outside radius in meters
radiusIn =  radiusOut * exp(-Z0 / 60);   % inside radius on meters
voltageIn  = 1;                          % voltage of inside conductor
voltageOut = 0;                          % voltage of outside conductor

refine = 2;                              % create denser mesh packed with triangles
jiggle = 'off';                          % jiggle mesh to spread nodes more evenly

[points, edges, triangles, geom] = geomcable(radiusIn, radiusOut, refine, jiggle);

numNodes = size(points, 2);              % total number of nodes in mesh
numBoarder = size(edges, 2);             % total number of edges at the boarder of mesh
numElements = size(triangles, 2);        % total number of triangle elements



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Boundary Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


InitValues = zeros(numNodes,1);                  % contains the Initial Potential Values at each node 
nodesKnown = zeros(numNodes,1);                  % contains 1 if the node has known value, else 0
nodesBoundary = union(edges(1, :), edges(2, :)); % take indexes of nodes that lie at the boundaries

radius = vecnorm(points(:, nodesBoundary));      % radius of each node at the boundary
meanRadius = mean([radiusOut radiusIn]);         % find the mean of radiusIn|Out
isOut = nodesBoundary(radius > meanRadius);      % take the nodes at the outter boundary
isIn =  nodesBoundary(radius < meanRadius);      % take the nodes at the inner boundary

nodesKnown(nodesBoundary) = 1;                   % nodes at the boundaries have known values
InitValues(isOut) = voltageOut;                  % assign Dirichlet conditions at the outher boundary
InitValues(isIn) =  voltageIn;                   % assign Dirichlet conditions at the inner boundary
numUknowns = numNodes - sum(nodesKnown);         % degrees of freedom



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PDE solver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


method = 'direct';                          % method of solver direct, bicg, gmres
e0 = 8.854 * 1e-12;                         % vacuum permittivity
permittivity = e0 * ones(1,numElements);    % permittivity inside capacitor

[Potential, StiffnessGlobal] = pdesolver(points, triangles, permittivity, nodesKnown, InitValues, method);


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


VoltDiff = voltageIn - voltageOut;                            % voltage difference in the 2 conductors
Energy = 1/2 * Potential' * (StiffnessGlobal * Potential);    % the approximated energy of the capacitor
CapacitancePredict = 2 * Energy / VoltDiff ^ 2;               % the approximated capacitance
CapacitanceTheory =  2 * pi * e0 / log(radiusOut / radiusIn); % analytic calculation of capacitance

error = abs(CapacitanceTheory - CapacitancePredict);          % error of aproximation
percent = 100 * error / CapacitanceTheory;                    % percentege of error in relation to analytic capacitance

