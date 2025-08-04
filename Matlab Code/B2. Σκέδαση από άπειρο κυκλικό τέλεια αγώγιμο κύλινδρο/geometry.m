function [points, edges, triangles, geom] = ...
geometry(scatRadius, pmlDepth, airGap, refine, jiggle)
%************************************************************************
% This function generates the triangle mesh for a specified geometry
%
% Input Arguments:
% scatter_radius radius of the infinite conductive cylinder
% pml_depth      thickness of pml layer
% air_gap        gap between pml and scatter
% refine         number of mesh refinements
% jiggle         jiggle the mesh to improve quality
%
% Output Arguments:
% points         the mesh point nodes
% edges          the mesh edges
% triangles      the mesh triangles
% geom           decomposed solid geometry matrix
%************************************************************************
for i = 1 : 1
    if scatRadius <= 0
        error('Error. \nRadius %s=%.3f must be greater than 0.', inputname(1),radius);
    
    elseif pmlDepth <= 0
        error('Error. \nPML Depth %s=%.3f must be greater than 0.', inputname(2),pmlDepth);
        
    elseif airGap <= 0
        error('Error. \nAir gap %s=%.3f must be greater than 0.', inputname(3),airGap);
    
    elseif refine < 0
        refine = 0;
    
    elseif ~isequal(jiggle, 'on') && ~isequal(jiggle, 'off')
        error('Error. \nJiggle can be either on or off, not %s', jiggle);
    
    elseif nargin == 3
        
        % the defaults
        refine = 1;
        jiggle = 'on';
    end
end             % check for the correct format of input arguments



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry of the Scattering problem from an infitite cylindrical conductor
%                    with PML layers at the boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


w = airGap + scatRadius;   
d = pmlDepth + w;                              
                                                             % subdomain numbering
scatterer =         [1 0 0 scatRadius 0 0 0 0 0 0]';     % F0   
compDomain =       [3 4 [-w w w -w] [-w -w w w]]';          % F1
pmlLeftX =         [2 4 [-d -w -w -d] [-w -w w w]]';        % F3
pmlRightX =        [2 4 [w d d w] [-w -w w w]]';            % F5
pmlUpY =           [2 4 [-w w w -w] [w w d d]]';            % F6
pmlDownY =         [2 4 [-w w w -w] [-d -d -w -w]]';        % F9
pmlLeftUpXY =      [2 4 [-d -w -w -d] [w w d d]]';          % F4
pmlLeftDownXY =    [2 4 [-d -w -w -d] [-d -d -w -w]]';      % F2
pmlRightUpXY =     [2 4 [w d d w] [w w d d]]';              % F7
pmlRightDownXY =   [2 4 [w d d w] [-d -d -w -w]]';          % F8


geom_descrition = [compDomain, scatterer, pmlLeftX, pmlRightX, pmlUpY, pmlDownY,...
                   pmlLeftUpXY, pmlLeftDownXY, pmlRightUpXY, pmlRightDownXY];           

nameSpace = char('CD','SC','pL','pR','pU','pD','pLU','pLD','pRU','pRD')';                       
setFormula = 'CD - SC + pL + pR + pU + pD + pLU + pLD + pRU + pRD';   

geom = decsg(geom_descrition, setFormula, nameSpace);    % decomposed solid geometry matrix
[points, edges, triangles] = initmesh(geom);               % create the triangle mesh



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Create additional refinements and appling jiggle mesh at every
%        iteration if wanted. Also you specify in which region to apply 
%        the refinement using the Face number at the end of refinemesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1 : refine   
    
    [points, edges, triangles] = refinemesh(geom, points, edges, triangles); 
    
    % jiggle the mesh   
    if isequal(jiggle, 'on')         
        points = jigglemesh(points,edges,triangles,'opt','mean','iter',inf);
    end
end

