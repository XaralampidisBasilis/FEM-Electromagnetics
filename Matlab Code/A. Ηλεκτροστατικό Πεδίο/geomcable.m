function [points, edges, triangles, geom] = geomcable(radiusIn, radiusOut,refine,jiggle)
%************************************************************************
% This function generates the mesh triplet for a coaxial cable capacitor
%
% Input Arguments:
% ra        radius of the inside cylinder
% rb        radius of the outside cylinder
% refine    number of mesh refinements
% jiggle    jiggle the mesh to improve quality
%
% Output Arguments:
% p         the mesh point nodes
% e         the mesh edges
% tr        the mesh triangles
%************************************************************************
for i = 1 : 1
    if radiusIn >= radiusOut
        error('Error. \nRadius %s=%.3f must be greater than %s=%.3f.', inputname(2),radiusOut,inputname(1),radiusIn);

    elseif refine < 0
        refine = 0;
    
    elseif ~isequal(jiggle, 'on') && ~isequal(jiggle, 'off')
        error('Error. \nJiggle can be either on or off, not %s', jiggle);
    end
end         % check for the correct format of input arguments


circleA = [1 0 0 radiusIn]';   % circle of center (0,0) and radius a
circleB = [1 0 0 radiusOut]';   % circle of center (0,0) and radius b
gd = [circleA circleB];  % geometry description

setFormula = 'R2 - R1';          % exclude the inside circle from the outer
nameSpace = char('R1','R2')';   % name the regions R1 R2
geom = decsg(gd,setFormula,nameSpace);    % decomposed solid geometry matrix

[points,edges,triangles] = initmesh(geom);  % create the triangle mesh


% do any number of refinements
for i = 1 : refine   
    [points,edges,triangles] = refinemesh(geom,points,edges,triangles,1); 
    % jiggle the mesh   
    if isequal(jiggle, 'on')         
        points = jigglemesh(points,edges,triangles,'opt','mean','iter',inf);
    end
end
