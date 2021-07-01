function [BCMatrix, BCRHS] = boundaryCondition3D(BC)
% It creates the matrix of coefficient based on the BC structureprovided
% by the user. It also generates the right hand side vector of the linear
% system of equations
%
% SYNOPSIS:
%   [BCMatrix, BCRHS] = boundaryCondition3D(BC)
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

% Note: I use a for loop here for more readability of the code!

% extract data from the mesh structure
Nxyz = BC.domain.dims;
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
dx_1 = BC.domain.cellsize.x(1);
dx_end = BC.domain.cellsize.x(end);
dy_1 = BC.domain.cellsize.y(1);
dy_end = BC.domain.cellsize.y(end);
dz_1 = BC.domain.cellsize.z(1);
dz_end = BC.domain.cellsize.z(end);


% number of boundary nodes (axact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 8*((Nx+1)*(Ny+1)+(Nx+1)*(Nz+1)+(Ny+1)*(Nz+1));

% define the vectors to be used for the creation of the sparse matrix
ii = zeros(nb,1);
jj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
BCRHS = zeros((Nx+2)*(Ny+2)*(Nz+2), 1);

% assign value to the corner nodes (useless cells)
q = 1:8;
ii(q) = BC.domain.corners; jj(q) = BC.domain.corners;
s(q) = 1; BCRHS(BC.domain.corners) = 0;

% assign values to the edges (useless cells)
q = q(end)+(1:length(BC.domain.edges));
ii(q) = BC.domain.edges; jj(q) = BC.domain.edges;
s(q) = 1; BCRHS(BC.domain.edges) = 0;

% Assign values to the boundary condition matrix and the RHS vector based
% on the BC structure
if (BC.top.periodic ==0) && (BC.bottom.periodic == 0)
    % top boundary
    j=Ny+2;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.top.b/2 + BC.top.a/dy_end;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j-1,k); s(q) = BC.top.b/2 - BC.top.a/dy_end;
    BCRHS(G(i,j,k)) = BC.top.c;

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j+1,k);  s(q) = -(BC.bottom.b/2 + BC.bottom.a/dy_1);
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.bottom.b/2 - BC.bottom.a/dy_1);
    BCRHS(G(i,j,k)) = -(BC.bottom.c);
elseif (BC.top.periodic ==1) || (BC.bottom.periodic == 1) % periodic
    % top boundary
    j=Ny+2;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j-1,k);  s(q) = -1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,1,k); s(q) = dy_end/dy_1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,2,k); s(q) = -dy_end/dy_1;
    BCRHS(G(i,j,k)) = 0;

    % Bottom boundary
    j=1;
    i=2:Nx+1;
    k=2:Nz+1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j+1,k);  s(q) = 1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,Ny+1,k); s(q) = -1;
    q = q(end)+(1:Nx*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,Ny+2,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

if (BC.right.periodic == 0) && (BC.left.periodic == 0)
    % Right boundary
    i=Nx+2;
    j=2:Ny+1;
    k=2:Nz+1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.right.b/2 + BC.right.a/dx_end;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i-1,j,k); s(q) = BC.right.b/2 - BC.right.a/dx_end;
    BCRHS(G(i,j,k)) = BC.right.c;

    % Left boundary
    i = 1;
    j=2:Ny+1;
    k=2:Nz+1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i+1,j,k);  s(q) = -(BC.left.b/2 + BC.left.a/dx_1);
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.left.b/2 - BC.left.a/dx_1);
    BCRHS(G(i,j,k)) = -(BC.left.c);
elseif (BC.right.periodic == 1) || (BC.left.periodic == 1) % periodic
    % Right boundary
    i=Nx+2;
    j=2:Ny+1;
    k=2:Nz+1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i-1,j,k);  s(q) = -1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(1,j,k); s(q) = dx_end/dx_1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(2,j,k); s(q) = -dx_end/dx_1;
    BCRHS(G(i,j,k)) = 0;

    % Left boundary
    i = 1;
    j=2:Ny+1;
    k=2:Nz+1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(i+1,j,k);  s(q) = 1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(Nx+1,j,k); s(q) = -1;
    q = q(end)+(1:Ny*Nz);
    ii(q) = G(i,j,k);  jj(q) = G(Nx+2,j,k); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

if (BC.front.periodic == 0) && (BC.back.periodic == 0)
    % Front boundary
    k=Nz+2;
    i = 2:Nx+1;
    j=2:Ny+1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = BC.front.b/2 + BC.front.a/dz_end;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k-1); s(q) = BC.front.b/2 - BC.front.a/dz_end;
    BCRHS(G(i,j,k)) = BC.front.c;

    % Back boundary
    k=1;
    i = 2:Nx+1;
    j=2:Ny+1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k+1);  s(q) = -(BC.back.b/2 + BC.back.a/dz_1);
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k); s(q) = -(BC.back.b/2 - BC.back.a/dz_1);
    BCRHS(G(i,j,k)) = -(BC.back.c);
elseif (BC.front.periodic == 1) || (BC.back.periodic == 1) % periodic
    % Front boundary
    k=Nz+2;
    i = 2:Nx+1;
    j=2:Ny+1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k-1);  s(q) = -1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,1); s(q) = dz_end/dz_1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,2); s(q) = -dz_end/dz_1;
    BCRHS(G(i,j,k)) = 0;

    % Back boundary
    k=1;
    i = 2:Nx+1;
    j=2:Ny+1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k);  s(q) = 1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,k+1);  s(q) = 1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,Nz+1); s(q) = -1;
    q = q(end)+(1:Nx*Ny);
    ii(q) = G(i,j,k);  jj(q) = G(i,j,Nz+2); s(q) = -1;
    BCRHS(G(i,j,k)) = 0;
end

% Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii(1:q(end)), jj(1:q(end)), s(1:q(end)), ...
    (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));
