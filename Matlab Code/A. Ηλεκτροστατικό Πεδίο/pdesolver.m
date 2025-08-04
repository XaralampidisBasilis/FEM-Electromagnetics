function [X,Sg] = pdesolver(p,tr,er,nodesKnown,X0,method)
%************************************************************************
% This function solves the electrostatic pde for a given mesh
%
% p         the mesh point positions of nodes
% tr        the mesh triangles
% er        the constant relative permittivity inside the capacitor
% node_id   binary information if a node has known value or not
% X0        information about the initial potential values
% method    method of pde solver direct, bicg, gmres
%
% Output Arguments:
% X         the potential value of every node
% Sg        the global stiffness matrix
%************************************************************************


numNodes = size(p,2);                    % total number of nodes in mesh
numElements = size(tr,2);                % total number of triangle elements
numKnown = sum(nodesKnown);              % total number of known nodes
numUnknown = numNodes - numKnown;        % total number of the unknown nodes


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


Sg = zeros(numNodes);          % global stiffness matrix
Sff = zeros(numUnknown);         % global unkown stiffness matrix 
Se = zeros(3, 3);        % local stiffness matrix
B = zeros(numUnknown,1);
[x, y, a, b, c, node] = deal(zeros(3,1));

for ie = 1 : numElements

    node(1:3) = tr(1:3,ie);         % indexes of nodes at the triangle ie
    x(1:3) = p(1,node(1:3));        % x coordinates of trinagle nodes
    y(1:3) = p(2,node(1:3));        % y coordinates of triangle nodes 
    

    De = det([ones(3,1), x, y]);    % area of the triangle element ie
    Ae = abs(De/2);     
    for n = 1 : 3                   % calculate the simplex coordinates variable values
        i = mod(n, 3) + 1;
        j = mod(n + 1, 3) + 1;
        a(n) = (x(i) * y(j) - x(j) * y(i)) / De;
        b(n) = (y(i) - y(j)) / De;
        c(n) = (x(j) - x(i)) / De;
    end
    
    % Initialize the global stiffness matrix 
    for i = 1 : 3
        for j = 1 : 3
            
            % calculate the local stiffness matrix 
            Se(i,j) = er(ie) * (b(i) * b(j) + c(i) * c(j)) * Ae;
            
            % calculate the global stiffness matrix
            Sg(node(i), node(j)) = Sg(node(i), node(j)) + Se(i,j);
            
            % calculate the uknown global stiffness matrix
            if nodesKnown(node(i)) == 0       % if node(i) is uknown
                if nodesKnown(node(j)) == 0   % if node(j) is uknown

                    Sff(newIndex(node(i)), newIndex(node(j))) = ...
                    Sff(newIndex(node(i)), newIndex(node(j))) + Se(i,j);
                else                       
                    B(newIndex(node(i))) = ...
                    B(newIndex(node(i))) - Se(i,j) * X0(node(j)); 
                end
            end
            
        end
    end
end


tol = 10^-5;
maxit = 10000;
if isequal(method, 'bicg')
    X = bicg(Sff,B,tol,maxit);     % Biconjugate gradient method
elseif isequal(method, 'gmres')
    X = gmres(Sff,B,5,tol,maxit);  % Generalized minimum residual method 
else
    X = Sff\B;                     % solve the pde system with direct solver
end

%{
index array maps all the uknown values of nodes into a new indexing ranging
from 1 : Nf. Tha array X has this type of indexing. Now inv array is the 
inverse function of index that maps values from 1 : Nf to the old indexing 
of nodes in the p array that holds also in X0.
%}
[newIndex, inv] = sort(newIndex);
inv(newIndex == 0) = []; 
X0(inv) = X(1 : numUnknown);               
X = X0;                 % now X holds all the values for the nodes 

end
