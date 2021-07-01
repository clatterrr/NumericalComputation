function [M, Mx, My, Mz] = diffusionTerm3D(D)
% This function uses the central difference scheme to discretize a 2D
% diffusion term in the form \grad . (D \grad \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
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

% extract data from the mesh structure
Nx = D.domain.dims(1);
Ny = D.domain.dims(2);
Nz = D.domain.dims(3);
G=reshape((1:(Nx+2)*(Ny+2)*(Nz+2)), Nx+2, Ny+2, Nz+2);
DX = repmat(D.domain.cellsize.x, 1, Ny, Nz);
DY = repmat(D.domain.cellsize.y', Nx, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = D.domain.cellsize.z;
DZ=repmat(DZ, Nx, Ny, 1);
dx = 0.5*(DX(1:end-1,:,:)+DX(2:end,:,:));
dy = 0.5*(DY(:,1:end-1,:)+DY(:,2:end,:));
dz = 0.5*(DZ(:,:,1:end-1)+DZ(:,:,2:end));

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
jjx = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
sx = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
iiy = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
jjy = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
sy = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
iiz = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
jjz = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
sz = zeros(3*(Nx+2)*(Ny+2)*(Nz+2),1);
mnx = Nx*Ny*Nz;	mny = Nx*Ny*Nz;   mnz = Nx*Ny*Nz;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
Dx = D.xvalue;
Dy = D.yvalue;
Dz = D.zvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
De = Dx(2:Nx+1,:,:)./(dx(2:Nx+1,:,:).*DX(2:Nx+1,:,:));
Dw = Dx(1:Nx,:,:)./(dx(1:Nx,:,:).*DX(2:Nx+1,:,:));
Dn = Dy(:,2:Ny+1,:)./(dy(:,2:Ny+1,:).*DY(:,2:Ny+1,:));
Ds = Dy(:,1:Ny,:)./(dy(:,1:Ny,:).*DY(:,2:Ny+1,:));
Df = Dz(:,:,2:Nz+1)./(dz(:,:,2:Nz+1).*DZ(:,:,2:Nz+1));
Db = Dz(:,:,1:Nz)./(dz(:,:,1:Nz).*DZ(:,:,2:Nz+1));

% calculate the coefficients for the internal cells
AE = reshape(De,mnx,1);
AW = reshape(Dw,mnx,1);
AN = reshape(Dn,mny,1);
AS = reshape(Ds,mny,1);
AF = reshape(Df,mnz,1);
AB = reshape(Db,mnz,1);
APx = reshape(-(De+Dw),mnx,1);
APy = reshape(-(Dn+Ds),mny,1);
APz = reshape(-(Df+Db),mnz,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
rowz_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnz,1); % main diagonal z
iiz(1:3*mnz) = repmat(rowz_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nx,2:Ny+1,2:Nz+1),mnx,1); reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); reshape(G(3:Nx+2,2:Ny+1,2:Nz+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nx+1,1:Ny,2:Nz+1),mny,1); reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mny,1); reshape(G(2:Nx+1,3:Ny+2,2:Nz+1),mny,1)];
jjz(1:3*mnz) = [reshape(G(2:Nx+1,2:Ny+1,1:Nz),mnz,1); reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnz,1); reshape(G(2:Nx+1,2:Ny+1,3:Nz+2),mnz,1)];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sz(1:3*mnz) = [AB; APz; AF];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
kz = 3*mnz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));
M = Mx + My + Mz;
