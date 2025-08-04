function [points, edges, triangles, geom] = ...
geomplates(boxWidth, boxHeight, platesWidth, platesThickness, platesDistance, refine, jiggle)
%************************************************************************
% This function generates the mesh triplet for a parallel plates capacitor
% 
% Input Arguments:
% A         width of the computational box
% B         hight of the computetional box
% w         width of plates
% t         thicknes of plates
% d         distance between plates
% refine    number of mesh refinements
% jiggle    jiggle the mesh to improve quality
%
% Output Arguments:
% p         the mesh point nodes
% e         the mesh edges
% tr        the mesh triangles
%************************************************************************
for checkFormat = 1
    if refine < 0
        refine = 0;
    
    elseif ~isequal(jiggle, 'on') && ~isequal(jiggle, 'off')
        error('Error. \nJiggle can be either on or off, not %s', jiggle);
    end
end

A = boxWidth;
B = boxHeight;
w = platesWidth;
t = platesThickness;
d = platesDistance;
RA = platesWidth;
RB = platesWidth;

box =       [2 4 [-A/2 A/2 A/2 -A/2] [-B/2 -B/2 B/2 B/2]]';        % the computational box
plate1 =    [2 4 [-w/2 w/2 w/2 -w/2] [d/2 d/2 d/2+t d/2+t]]';      % the up plate
plate2 =    [2 4 [-w/2 w/2 w/2 -w/2] [-d/2 -d/2 -d/2-t -d/2-t]]';  % the down plate
between =   [2 4 [-w/2 w/2 w/2 -w/2] [-d/2 -d/2 d/2 d/2]]';        % the region between plates
refineBox = [2 4 [-RA RA RA -RA] [-RB -RB RB RB]]';                % the box in which we will aply the refinment

gd = [box refineBox, plate1 plate2, between];                  % geometry description
ns = char('Box','RefBox','PlateUp','PlateDown','Between')';    % name the regions R1 R2 R3 R4
sf = '(Box + RefBox + Between) - (PlateUp + PlateDown)';       % exclude the plates from the box


geom = decsg(gd,sf,ns);         % decomposed solid geometry matrix
[points,edges,triangles] = initmesh(geom);      % create the triangle mesh


% do any number of refinements
for i = 1 : refine   
            
    regions = [3 2];
    
    [points,edges,triangles] = refinemesh(geom,points,edges,triangles, regions);  
    
    % jiggle the mesh   
    if isequal(jiggle, 'on')         
        
        points = jigglemesh(points,edges,triangles,'opt','mean','iter',inf); 
    end
end

