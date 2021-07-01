function visualizeCellVectors3D(phi_cell)
%VISUALIZECELLS plots the values of cell variable phi
%
% SYNOPSIS:
%
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% Written by Ali A. Eftekhari
% See the license file

[X,Y,Z]=meshgrid(phi_cell.domain.cellcenters.y, phi_cell.domain.cellcenters.x, ...
    phi_cell.domain.cellcenters.z);
quiver3(X,Y,Z,phi_cell.xvalue, phi_cell.yvalue, phi_cell.zvalue);
xlabel('Cell centers [y vlaues]'); % this is correct [matrix not rotated]
ylabel('Cell centers [x vlaues]'); % this is correct [matrix not rotated]
zlabel('Cell centers [z vlaues]');
axis equal tight
% colorbar
