function pdeploter(points, edges, triangles, geom, Values, plotMesh)
%************************************************************************
%       This function plots the values of the geometry nodes
%************************************************************************


if isequal(plotMesh, 'on')
    
    % plot the created triangle mesh and geometry
    figure
    pdemesh(points, edges, triangles,'NodeLabels','off','ElementLabels','off');  
    hold on
    pdegplot(geom,'VertexLabels','off','EdgeLabels','off','FaceLabels','on');    
    hold off
    
end

% plot the TEM modes
figure
pdeplot(points,edges,triangles,'xydata',Values,...  % plot the node Values
        'colormap','jet');
colorbar('off');
axis tight;
axis equal;
axis off

end

