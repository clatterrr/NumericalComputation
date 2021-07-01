function M = convectionUpwindTermCylindrical1D(u)
% This function uses the upwind scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   M = convectionUpwindTermCylindrical1D(u)
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
Nr = u.domain.dims(1);
G = 1:Nr+2;
DXp = u.domain.cellsize.x(2:end-1);
rp = u.domain.cellcenters.x;
rf = u.domain.facecenters.x;

% define the vectors to store the sparse matrix data
iix = zeros(3*(Nr+2),1);
jjx = zeros(3*(Nr+2),1);
sx = zeros(3*(Nr+2),1);

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nr+1);		uw = ux(1:Nr);
re = rf(2:Nr+1);     rw = rf(1:Nr);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);

% calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(DXp.*rp),Nr,1);
AW = reshape(-rw.*uw_max./(DXp.*rp),Nr,1);
APx = reshape((re.*ue_max-rw.*uw_min)./(DXp.*rp),Nr,1);

% correct for the cells next to the boundary
% Left boundary:
APx(1) = APx(1)-rw(1)*uw_max(1)/(2*rp(1)*DXp(1));   AW(1) = AW(1)/2;
% Right boundary:
AE(end) = AE(end)/2;    APx(end) = APx(end)+re(end)*ue_min(end)/(2*rp(end)*DXp(end));

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1),Nr,1); % main diagonal x
iix(1:3*Nr) = repmat(rowx_index,3,1);
jjx(1:3*Nr) = [reshape(G(1:Nr),Nr,1); reshape(G(2:Nr+1),Nr,1); reshape(G(3:Nr+2),Nr,1)];
sx(1:3*Nr) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nr;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nr+2, Nr+2);
