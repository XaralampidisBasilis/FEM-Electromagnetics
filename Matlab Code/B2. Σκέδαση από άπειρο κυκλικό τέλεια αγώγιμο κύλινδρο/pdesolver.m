function EzScat = ...
pdesolver(points, triangles, ezz, mxx, myy, frequency, nodesKnown, InitValues)


omega = 2 * pi * frequency;                % angular frequency
numNodes = size(points,2);                 % total number of nodes in mesh
numElements = size(triangles,2);           % total number of triangle elements
numKnown = sum(nodesKnown);                % total number of known nodes
numUnknown =  numNodes - numKnown;         % total number of the unknown nodes


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


Se = zeros(3);                   % local stiffness matrix
Te = zeros(3);                   % local mass matrix
A = zeros(numUnknown);    
B = zeros(numUnknown, 1);

[x, y, a, b, c, nodesIdx] = deal(zeros(3,1));

for elementIdx = 1 : numElements

    nodesIdx(1:3) = triangles(1:3, elementIdx);     % indexes of nodes at the triangle ie
    x(1:3) = points(1, nodesIdx(1:3));                 % x coordinates of trinagle nodes
    y(1:3) = points(2, nodesIdx(1:3));                 % y coordinates of triangle nodes 
    

    De = det([ones(3,1), x, y]);    % area of the triangle element
    Ae = abs(De/2);     
    
    for n = 1 : 3                   % calculate the simplex coordinates variable values
        i = mod(n, 3) + 1;
        j = mod(n + 1, 3) + 1;
        a(n) = (x(i) * y(j) - x(j) * y(i)) / De;
        b(n) = (y(i) - y(j)) / De;
        c(n) = (x(j) - x(i)) / De;
    end
    
    for i = 1 : 3
        for j = 1 : 3
            
            % calculate the local matrices
            Se(i,j) = Ae *(1 / myy(elementIdx) * b(i) * b(j) + ...
                           1 / mxx(elementIdx) * c(i) * c(j));
            Te(i,j) = Ae * ezz(elementIdx) / 12;
            if i == j, Te(i,i) = 2 * Te(i,j);
            end
            
            
            % calculate the uknown matricies
            if nodesKnown(nodesIdx(i)) == 0       % if node(i) is uknown
                if nodesKnown(nodesIdx(j)) == 0   % if node(j) is uknown
                
                    A(newIndex(nodesIdx(i)), newIndex(nodesIdx(j))) = ...
                    A(newIndex(nodesIdx(i)), newIndex(nodesIdx(j))) + ...
                    Se(i,j) - omega^2 * Te(i,j);

                else
                    B(newIndex(nodesIdx(i))) = ...
                    B(newIndex(nodesIdx(i))) - ...
                    (Se(i,j) - omega^2 * Te(i,j)) * InitValues(nodesIdx(j));

                end
     
            end
            
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                PDE solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EzScatff = A \ B;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Gather all the values of the nodes and return to the old indexing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[newIndex, oldIndex] = sort(newIndex);
oldIndex(newIndex == 0) = []; 
InitValues(oldIndex) = EzScatff(1 : numUnknown);  
EzScat = InitValues;                     

end
