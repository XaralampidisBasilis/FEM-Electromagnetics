function [points, edges, triangles, geom] = geometry(radius, refine, jiggle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the mesh triplet for a cylindrical hollow waveguide 
%
% Input:
%
%   radius     radius of the cylindical waveguide
%   refine     number of mesh refinements
%   jiggle     jiggle the mesh to improve quality
%
% Output:
%
%   points     coordinates of the mesh point nodes
%   edges      matrix desctibing the mesh edges
%   triangles  matrix describing the mesh triangles
%   geom       decomposed solid geometry matrix
%
% Format:
% 
%   radius =    positive scalar number
%   refine =    natural number
%   jiggle =    either 'on' of 'off'
%   points =    matrix of size [2 n] for n-points
%   eges   =    matrix of size [7 m] for m-edges
%   triangles = matrix of size [4 e] for e-triangles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for checkFormat = 1
    
    if radius <= 0
        error('Error. Radius must be > 0, not %s',radius);

    elseif refine < 0
        disp('refine became zero');
        refine = 0;
    
    elseif ~isequal(jiggle, 'on') && ~isequal(jiggle, 'off')
        error('Error. Jiggle can be either on or off, not %s', jiggle);
    end
end         


circle = [1 0 0 radius 0 0 0 0 0 0];   % circle of center (0,0) and specified radius 
box1    = [2 4 [-radius radius radius -radius] [0 0 radius radius]];       % added only to create a line at the middle
box2    = [2 4 [-radius radius radius -radius] [-radius -radius 0 0]];     % it is not needed, just wanted to break the symetry

geomDescrition = [circle', box1', box2'];    % geometry description

nameSpace = char('C','B1','B2')';                     % name the existing regions
setFormula = '(C * B1) + (C * B2)';                   % describe how the regions will coexist
geom = decsg(geomDescrition, setFormula, nameSpace);  % decomposed solid geometry matrix

[points, edges, triangles] = initmesh(geom);          % create the triangle mesh


% do any number of refinements
for iterations = 1 : refine   
    
    [points, edges, triangles] = refinemesh(geom, points, edges, triangles); 
    
    % jiggle the mesh   
    if isequal(jiggle, 'on')         
        points = jigglemesh(points,edges,triangles,'opt','mean','iter',inf);
    end
end


