clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 B.1 Cylindrical Dielectric Waveguide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


radius = 0.01;                     % radius of the waveguide in meters   
numModes = 9;                     % the number of unique modes we want to get 
modeType = 'TM';                   % TM(Hz = 0) or TE(Ez = 0) solutions of the eigenvalue problem

refine = 3;                        % create denser mesh packed with triangles
jiggle = 'off';                    % jiggle mesh to spread nodes more evenly

[points, edges, triangles, geom] = geometry(radius, refine, jiggle);

numNodes =    size(points, 2);            % total number of nodes in mesh
numBoarder =  size(edges, 2);             % total number of edges at the boarder of mesh
numElements = size(triangles, 2);         % total number of triangle elements


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Boundary Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


InitValues = zeros(numNodes, 1);          % the initial values of the nodes
nodesKnown = zeros(numNodes, 1);          % contains 1 if the node has known value, else 0
numUknowns = numNodes - sum(nodesKnown);  % degrees of freedom

if isequal(modeType, 'TM')
    
    edgeBoundary = pdesde(edges, [1 2]);      % take indexes of edges that lie at
                                              % the boudary of subdomain 1 
    nodesBoundary = edges(1, edgeBoundary);   % take the nodes at this boundary
    nodesKnown(nodesBoundary) = 1;            % nodes at the boundary have known values
    InitValues(nodesBoundary) = 0;            % assign Dirichlet conditions for TM modes
                                              % and Neuman for TE modes
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PDE solver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numEigenVals = round(2.5 * numModes);     % the number of modes we will get from solving 
                                        % the eigenvalue problem. But someome modes are 
                                        % not unique and they repeat so we need to filter 
                                        % those out

[eigVec, eigVal] = ...
pdesolver(points, triangles, nodesKnown, InitValues, numEigenVals);

                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Search for artificial non unique modes and remove them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for TE modes we get a zero eigen value 
% which is not an acceptable solution 
if isequal(modeType, 'TE') && abs(eigVal(1)) <= 1e-5
    eigVal(1) = [];
    eigVec(:, 1) = [];
end

[error, percent, Fields, cutofFreqTheory] = ... 
errorfinder(eigVec, eigVal, numModes, modeType, radius);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PDE ploter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotMesh = 'off';
plotResults = 'on';

if isequal(plotResults, 'on')
    for mode = 1 : numModes 

        X = Fields(:, mode);
        pdeploter(points, edges, triangles, geom, X, plotMesh);
        saveas(gcf,[modeType,'n', num2str(mode),'r',num2str(refine),'.png']);
        close(gcf);

    end
end
