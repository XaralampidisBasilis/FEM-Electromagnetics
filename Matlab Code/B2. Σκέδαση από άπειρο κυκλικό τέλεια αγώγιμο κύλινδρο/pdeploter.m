function pdeploter(points, edges, triangles, geom, values)
%************************************************************************
%       This function plots the values of the geometry nodes
%************************************************************************

% plot the created triangle mesh and geometry
figure
pdemesh(points, edges, triangles,'NodeLabels','off','ElementLabels','off');  
hold on
pdegplot(geom,'VertexLabels','off','EdgeLabels','off','FaceLabels','on');    
hold off


% get only the nodes inside the computational domain
nodes = 1 : size(points, 2);
[nodesIn, nodesBoundary] = pdesdp(points, edges, triangles, 1);
nodesCD = union(nodesIn, nodesBoundary);

nodesPML = setdiff(nodes, nodesCD);
values(nodesPML) = 0;


% plot the Amplitute of the field inside the computational domain
figure
pdeplot(points,edges,triangles,'xydata',values, ...  
        'colormap','jet');
axis tight;
axis equal;

end

