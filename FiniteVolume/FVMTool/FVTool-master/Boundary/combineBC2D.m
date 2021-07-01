function [Mout, RHSout] = combineBC2D(BC, Meq, RHSeq)
%COMBINEBC This function combines the boundary condition equations with the
%main physical model equations, and delivers the matrix of coefficient and
%RHS to be solved for the internal cells.
%
% SYNOPSIS:
%   [Mout, RHSout] = combineBC2D(BC, Meq, RHSeq)
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

% extract data from the mesh structure
Nxy = BC.domain.dims;
Nx = Nxy(1); Ny = Nxy(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
dx = BC.domain.cellsize.x;
dy = BC.domain.cellsize.y;

% define the RHS column vector
ms = size(Meq);
M = Meq;
RHS = RHSeq;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
% top boundary
j=Ny+2;
i=2:Nx+1;
top = reshape(sub2ind(ms, G(i,j-1), G(i,j-1)), Nx,1); % top boundary cells
topN = reshape(sub2ind(ms, G(i,j-1), G(i,j)), Nx, 1); % north cells to top boundary cells
M(top) = M(top)-((BC.top.b/2 - BC.top.a/dy(end))./(BC.top.b/2 + BC.top.a/dy(end))).*M(topN);
RHS(G(i,j-1)) = RHS(G(i,j-1))-M(topN).*BC.top.c./(BC.top.b/2 + BC.top.a/dy(end));

% Bottom boundary
j=1;
i=2:Nx+1;
bottom = reshape(sub2ind(ms, G(i,j+1), G(i,j+1)), Nx,1); % bottom boundary cells
bottomS = reshape(sub2ind(ms, G(i,j+1), G(i,j)), Nx, 1); % south cells to bottom boundary cells
M(bottom) = M(bottom)-((BC.bottom.b/2 + BC.bottom.a/dy(1))./(BC.bottom.b/2 - BC.bottom.a/dy(1))).*M(bottomS);
RHS(G(i,j+1)) = RHS(G(i,j+1))-M(bottomS).*BC.bottom.c./(BC.bottom.b/2 - BC.bottom.a/dy(1));

% Right boundary
i=Nx+2;
j=2:Ny+1;
right = reshape(sub2ind(ms, G(i-1,j), G(i-1,j)), Ny,1); % right boundary cells
rightE = reshape(sub2ind(ms, G(i-1,j), G(i,j)), Ny, 1); % east cells to right boundary cells
M(right) = M(right)-((BC.right.b/2 - BC.right.a/dx(end))./(BC.right.b/2 + BC.right.a/dx(end)))'.*M(rightE);
RHS(G(i-1,j)) = RHS(G(i-1,j))-M(rightE).*(BC.right.c./(BC.right.b/2 + BC.right.a/dx(end)))';

% Left boundary
i = 1;
j=2:Ny+1;
left = reshape(sub2ind(ms, G(i+1,j), G(i+1,j)), Ny,1); % left boundary cells
leftW = reshape(sub2ind(ms, G(i+1,j), G(i,j)), Ny, 1); % west cells to left boundary cells
M(left) = M(left)-((BC.left.b/2 + BC.left.a/dx(1))./(BC.left.b/2 - BC.left.a/dx(1)))'.*M(leftW);
RHS(G(i+1,j)) = RHS(G(i+1,j))-M(leftW).*(BC.left.c./(BC.left.b/2 - BC.left.a/dx(1)))';

Mout = M(G(2:end-1,2:end-1), G(2:end-1,2:end-1));
RHSout = RHS(reshape(G(2:end-1,2:end-1),Nx*Ny,1));
