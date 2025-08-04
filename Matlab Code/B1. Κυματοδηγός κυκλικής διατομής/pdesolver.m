function [eigVec, eigVal] = ...
pdesolver(points, triangles, nodesKnown, InitValues, numEigenVals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves the eigenvalue problem for the hollow cylindrical wavegiude 
%
% Input:
%   points         The coordinates of mesh nodes
%   triangles      The mesh triangles
%   permitivity    Permittivity inside the cylindrical waveguide
%   nodesKnown     Binary information if a node has known value 1 or not 0
%   InitValues     Information about the initial values
%   numEigenVals   The number of solutions we want
%
% Output:
%
%   eigVec         the values of the fields at every node
%   eigVal         frequencies at which electromagnetic modes appear
%   error          the error of cutofFreq from the theoretical predictions
%   percent        the same error but in percentage of the theoretical 
%
% Format:
% 
%   points =        matrix of size [2 n] for n-points
%   triangles =     matrix of size [4 e] for e-triangles
%   permitivity =   matrix of size [4 e] for e-triangles
%   nodesKnown =    matrix of size [n 1] for n-points
%   InitValues =    matrix of size [n 1] for n-points
%   numEigenVals =  natural number
%   modeType =      either 'TM' of 'TE'
%   eigVec =        matrix of size [n, maxModes] for n-points
%   eigVal =        matrix of size [maxModes 1] 
%   error =         matrix of size [maxModes 1]
%   percent =       matrix of size [maxModes 1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numNodes =    size(points,2);                % total number of nodes in mesh
numElements = size(triangles,2);             % total number of triangle elements
numUknown =   numNodes - sum(nodesKnown);    % total number of the unknown nodes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Assign new indexing to all the unknown nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter = 0;                  
newIndex = zeros(numNodes,1);    

for oldIndex = 1 : numNodes
    if nodesKnown(oldIndex) == 0        
        counter = counter + 1; 
        newIndex(oldIndex) = counter;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the discrete matricies for the Galerkin method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Se = zeros(3);                  % local stiffness matrix
Te = zeros(3);                  % local mass matrix
Sff = zeros(numUknown);         % uknown stiffness matrix
Tff = zeros(numUknown);         % uknown mass matrix
nodesIdx = zeros(3, 1);         % indicies of nodes at triangle element
X = zeros(3, 1);                % x coordinates of trinagle nodes
Y = zeros(3, 1);                % y coordinates of triangle nodes 
[a, b, c] = deal(zeros(3,1));   % simplex coordinates

for elementIdx = 1 : numElements

    nodesIdx(1:3) = triangles(1:3, elementIdx);  
    X(1:3) = points(1, nodesIdx(1:3));            
    Y(1:3) = points(2, nodesIdx(1:3));             
    De = det([ones(3,1), X, Y]);                            
    Ae = abs(De/2);              % area of the triangle element
    

    % calculate the simplex coordinates variable values
    for n = 1 : 3                   
        i = mod(n, 3) + 1;
        j = mod(n + 1, 3) + 1;
        a(n) = (X(i) * Y(j) - X(j) * Y(i)) / De;
        b(n) = (Y(i) - Y(j)) / De;
        c(n) = (X(j) - X(i)) / De;
    end
    
    for i = 1 : 3
        for j = 1 : 3
            
            % calculate the local matrices
            Se(i,j) = (b(i) * b(j) + c(i) * c(j)) * Ae;
            Te(i,j) = Ae / 12;
            if i == j, Te(i,i) = 2 * Te(i,j);
            end
            
            % calculate the uknown matricies
            if nodesKnown(nodesIdx(i)) == 0       % if node(i) is uknown
                if nodesKnown(nodesIdx(j)) == 0   % if node(j) is uknown

                    Sff(newIndex(nodesIdx(i)), newIndex(nodesIdx(j))) = ...
                    Sff(newIndex(nodesIdx(i)), newIndex(nodesIdx(j))) + Se(i,j);
                
                    Tff(newIndex(nodesIdx(i)), newIndex(nodesIdx(j))) = ...
                    Tff(newIndex(nodesIdx(i)), newIndex(nodesIdx(j))) + Te(i,j);
                end
            end
            
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the eigenvalue problem, taking the first k smallest eigenvalues and
% eigen vectors. V has the eigenvectors in columns and D is diagonal with
% the eigenvalues. Also want to make sure to keep only the unique modes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply small pertubation to get unique solutions at Neuman boundary problem
% https://scicomp.stackexchange.com/questions/25449/integrate-result-of-finite-element-calculation-in-matlab
pertubation = 1e-15;

[eigVec, eigVal] = eigs(Sff - pertubation, Tff, numEigenVals, 'sm');
eigVal = abs(real(diag(eigVal)));
eigVec = real(eigVec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Gather all the values of the nodes and return to the old indexing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[newIndex, oldIndex] = sort(newIndex);
oldIndex(newIndex == 0) = []; 
InitValues(oldIndex, 1 : numEigenVals) = eigVec(1 : numUknown, 1 : numEigenVals);  
eigVec = InitValues;                     


end
