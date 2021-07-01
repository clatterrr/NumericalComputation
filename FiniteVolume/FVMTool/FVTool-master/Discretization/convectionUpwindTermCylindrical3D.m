function [M, Mx, My, Mz] = convectionUpwindTermCylindrical3D(u)
% This function uses the upwind scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
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
Nr = u.domain.dims(1);
Ntetta = u.domain.dims(2);
Nz = u.domain.dims(3);
G=reshape((1:(Nr+2)*(Ntetta+2)*(Nz+2)), Nr+2, Ntetta+2, Nz+2);
DRp = repmat(u.domain.cellsize.x(2:end-1), 1, Ntetta, Nz);
DTHETAp = repmat(u.domain.cellsize.y(2:end-1)', Nr, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = u.domain.cellsize.z;
DZp=repmat(DZ(1,1,2:end-1), Nr, Ntetta, 1);
rp = repmat(u.domain.cellcenters.x, 1, Ntetta, Nz);
rf = repmat(u.domain.facecenters.x, 1, Ntetta, Nz);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
jjx = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
sx = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
iiy = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
jjy = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
sy = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
iiz = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
jjz = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
sz = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
mnx = Nr*Ntetta*Nz;	mny = Nr*Ntetta*Nz;   mnz = Nr*Ntetta*Nz;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;
uz = u.zvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nr+1,:,:);		uw = ux(1:Nr,:,:);
vn = uy(:,2:Ntetta+1,:);     vs = uy(:,1:Ntetta,:);
wf = uz(:,:,2:Nz+1);     wb = uz(:,:,1:Nz);
re = rf(2:Nr+1,:,:);         rw = rf(1:Nr,:,:);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);
vn_min = min(vn,0);	vn_max = max(vn,0);
vs_min = min(vs,0);	vs_max = max(vs,0);
wf_min = min(wf,0);	wf_max = max(wf,0);
wb_min = min(wb,0);	wb_max = max(wb,0);

% calculate the coefficients for the internal cells
AE = re.*ue_min./(DRp.*rp);
AW = -rw.*uw_max./(DRp.*rp);
AN = vn_min./(DTHETAp.*rp);
AS = -vs_max./(DTHETAp.*rp);
AF = wf_min./DZp;
AB = -wb_max./DZp;
APx = (re.*ue_max-rw.*uw_min)./(DRp.*rp);
APy = (vn_max-vs_min)./(DTHETAp.*rp);
APz = (wf_max-wb_min)./DZp;

% Also correct for the boundary cells (not the ghost cells)
% Left boundary:
APx(1,:,:) = APx(1,:,:)-rw(1,:,:).*uw_max(1,:,:)./(2*rp(1,:,:)*DRp(1));   AW(1,:,:) = AW(1,:,:)/2;
% Right boundary:
AE(end,:,:) = AE(end,:,:)/2;    APx(end,:,:) = APx(end,:,:)+re(end,:,:).*ue_min(end,:,:)./(2*DRp(end)*rp(end,:,:));
% Bottom boundary:
APy(:,1,:) = APy(:,1,:)-vs_max(:,1,:)./(2*DTHETAp(1)*rp(:,1,:));   AS(:,1,:) = AS(:,1,:)/2;
% Top boundary:
AN(:,end,:) = AN(:,end,:)/2;    APy(:,end,:) = APy(:,end,:)+vn_min(:,end,:)./(2*DTHETAp(end)*rp(:,end,:));
% Back boundary:
APz(:,:,1) = APz(:,:,1)-wb_max(:,:,1)/(2*DZp(1));   AB(:,:,1) = AB(:,:,1)/2;
% Front boundary:
AF(:,:,end) = AF(:,:,end)/2;    APz(:,:,end) = APz(:,:,end) + wf_min(:,:,end)/(2*DZp(end));

AE = reshape(AE,mnx,1);
AW = reshape(AW,mnx,1);
AN = reshape(AN,mny,1);
AS = reshape(AS,mny,1);
AF = reshape(AF,mnz,1);
AB = reshape(AB,mnz,1);
APx = reshape(APx,mnx,1);
APy = reshape(APy,mny,1);
APz = reshape(APz,mnz,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
rowz_index = reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnz,1); % main diagonal z
iiz(1:3*mnz) = repmat(rowz_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nr,2:Ntetta+1,2:Nz+1),mnx,1); reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnx,1); reshape(G(3:Nr+2,2:Ntetta+1,2:Nz+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nr+1,1:Ntetta,2:Nz+1),mny,1); reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mny,1); reshape(G(2:Nr+1,3:Ntetta+2,2:Nz+1),mny,1)];
jjz(1:3*mnz) = [reshape(G(2:Nr+1,2:Ntetta+1,1:Nz),mnz,1); reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnz,1); reshape(G(2:Nr+1,2:Ntetta+1,3:Nz+2),mnz,1)];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sz(1:3*mnz) = [AB; APz; AF];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
kz = 3*mnz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nr+2)*(Ntetta+2)*(Nz+2), (Nr+2)*(Ntetta+2)*(Nz+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nr+2)*(Ntetta+2)*(Nz+2), (Nr+2)*(Ntetta+2)*(Nz+2));
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), (Nr+2)*(Ntetta+2)*(Nz+2), (Nr+2)*(Ntetta+2)*(Nz+2));
M = Mx + My + Mz;