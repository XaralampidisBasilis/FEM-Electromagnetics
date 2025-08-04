clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     B.2 Electromagnetic Scattering from an Infinite Conductive Cylinder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


freq = 300 * 1e+9;                            % frequency of the EM wave
wavelength = physconst('LightSpeed') / freq;  % wavelength of the EM wave in air
scatRadius = wavelength * [1/2, 2, 5];        % radius of the infinite conductive cylinder
pmlDepth = wavelength * [1, 1/4];             % thickness of pml layer
airGap = wavelength;                          % gap between pml and scatter

refine = 3;                                % create denser mesh packed with triangles
jiggle = 'off';                            % jiggle mesh to spread nodes more evenly

scatRadius = scatRadius(1);
pmlDepth = pmlDepth(1);

[points, edges, triangles, geom] = geometry(scatRadius, pmlDepth, airGap, refine, jiggle);

numNodes =    size(points, 2);            % total number of nodes in mesh
numBoarder =  size(edges, 2);             % total number of edges at the boarder of mesh
numElements = size(triangles, 2);         % total number of triangle elements



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Boundary Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


InitValues = zeros(numNodes, 1);          % the initial values of the nodes
nodesKnown = zeros(numNodes, 1);          % contains if a node has known value 1 or not 0
numUknowns =  numNodes - sum(nodesKnown); % degrees of freedom
edgeBoundary = pdesde(edges, 1);          % take indexes of edges that lie at
                                          % the boudary of subdomain 1 
nodesBoundary = edges(1, edgeBoundary);   % take the nodes at this boundary
nodesKnown(nodesBoundary) = 1;           

k0 = 2 * pi /wavelength;                     % wavenumber in the air
E0 = 1;  % V/m                               % amplitude of TEM incident wave
FzInc = @(x) E0 * exp(-1j * k0 * x);         % Equation for the incident TEM wave 
EzScat = -FzInc(points(1, nodesBoundary));   % On the surface of the scatterer
InitValues(nodesBoundary) = EzScat;          % we have non homogenous Dirichlet conditions
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PDE solver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


reflectCoeff = 1e-6;
[ezz, mxx, myy] = material(triangles, pmlDepth, wavelength, reflectCoeff);    

EzScat = ...
pdesolver(points, triangles, ezz, mxx, myy, freq, nodesKnown, InitValues);

                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             PDE ploter only inside the computational domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EzInc = FzInc(points(1,:)');
Ez = EzScat + EzInc;     
Az =  abs(Ez);            % amplitute of the total Ez 
phz = angle(Ez);          % phase of the total Ez
pdeploter(points, edges, triangles, geom, Az);

