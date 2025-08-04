function [ezz, mxx, myy] = material(triangles, pmlDepth, wavelength, reflectCoeff)


e0 = 8.854 * 1e-12;                   % vacuum permittivity
m0 = 4 * pi * 1e-7;                   % Vacuum permeability
k0 = 2 * pi / wavelength;             % wavenumber in the air

b = -log(reflectCoeff) / (2 * k0 * pmlDepth);  % loss factor 
a = 1 - 1j * b;

numElements = size(triangles, 2);    % total number of triangle elements
erzz = zeros(numElements,1);          
mrxx = zeros(numElements,1);    
mryy = zeros(numElements,1);

% relative permitivity inside computational domain
indexCD = pdesdt(triangles, 1);
erzz(indexCD) = 1;
mrxx(indexCD) = 1;
mryy(indexCD) = 1;

% inside the x pml region
indexPMLx = pdesdt(triangles, [3 5]);
erzz(indexPMLx) = a;
mrxx(indexPMLx) = 1/a;
mryy(indexPMLx) = a;

% inside the y pml region
indexPMLy = pdesdt(triangles, [6 9]);
erzz(indexPMLy) = a;
mrxx(indexPMLy) = 1/a;
mryy(indexPMLy) = a;

% inside the corners of pml
indexPMLxy = pdesdt(triangles, [4 2 7 8]);
erzz(indexPMLxy) = a * a;
mrxx(indexPMLxy) = 1;
mryy(indexPMLxy) = 1;

ezz = e0 * erzz;  
mxx = m0 * mrxx;
myy = m0 * mryy;

end
