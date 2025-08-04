function pdeploter(points, edges, triangles, geom, Potential, plotMesh, plotField, plotPotential)
%************************************************************************
% This function plots 3 figures. The mesh, the electic vector field and 
% the potential of the capacitor
%
% Input Arguments:
% p         the mesh point positions of nodes
% e         the edges of the mesh
% tr        the mesh triangles
% dl        decomposed geometry matrix
% X         the potential values of capacitor
%************************************************************************


[Ex,Ey] = pdegrad(points, triangles, Potential);     % calculate the gradient 
ElectricField = [Ex; Ey];                            % Electric vector field
Amplitute = vecnorm(ElectricField);                  % the Amplitute of the vector field

if isequal(plotMesh, 'on')
    
    figure
    pdemesh(points, edges, triangles,'NodeLabels','off','ElementLabels','off');  
    hold on
    pdegplot(geom,'VertexLabels','off','EdgeLabels','off','FaceLabels','on'); % plot the geometry of the capacitor   
    hold off
    axis tight;
    axis equal;
    
end
if isequal(plotField, 'on')
    
    % plot the Electric Amplitude
    f = figure;
    f.Color = '[1 1 1]';
    Value = log(Amplitute+1).^4.5;                  % just for visual normalize
    Value = max(Amplitute) * Value ./ max(Value);   % get the same scale of numbers
    s1 = pdeplot(points,edges,triangles,'xydata',Value,'Contour','off');  
    alpha(s1, 1);
    
    
    % plot the geometry of the capacitor   
    
    hold on
    s2 = pdegplot(geom,'VertexLabels','off','EdgeLabels','off','FaceLabels','off'); 
    s2.Color = [0.2, 0.2, 0.2];
    hold off
    
           
    % plot the Flectric vector field 
    hold on
    ElectricField = ElectricField ./ Amplitute;                    % normalize the vector field 
    ElectricField = (Amplitute.^0.6) .* ElectricField;             % make the vector field look better
    q = pdeplot(points,edges,triangles,'flowdata',ElectricField);  
    q.Color = 'w';
    q.LineWidth = 0.6;
    hold off
    %}
    axis tight;
    axis equal;
    colormap bone;  % jet, parula, pink, bone
    set(gca,'Color','[1 1 1');

end
if isequal(plotPotential, 'on')
    
    figure
    pdeplot(points,edges,triangles,'xydata',Potential,...  % plot the potential
                   'colormap','parula');
    axis tight;
    axis equal;
    
end
end

