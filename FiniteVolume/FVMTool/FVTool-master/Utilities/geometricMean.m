function phiFaceAverage = geometricMean(phi)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the geometric average on
% the cell faces, for a uniform mesh.
%
% SYNOPSIS:
%   phiFaceAverage = geometricMean(phi)
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

% extract data from the mesh structure

d = phi.domain.dimension;
if (d ==1) || (d==1.5) || (d==1.8)
    dx = phi.domain.cellsize.x;
    xvalue=exp((dx(1:end-1).*log(phi.value(1:end-1))+dx(2:end).*log(phi.value(2:end)))./(dx(2:end)+dx(1:end-1)));
    yvalue=[];
    zvalue=[];
elseif (d == 2) || (d == 2.5) || (d == 2.8)
    Nx = phi.domain.dims(1);
    Ny = phi.domain.dims(2);
    dx = repmat(phi.domain.cellsize.x, 1, Ny);
    dy = repmat(phi.domain.cellsize.y', Nx, 1);
	xvalue=exp((dx(1:end-1,:).*log(phi.value(1:end-1,2:end-1))+...
        dx(2:end,:).*log(phi.value(2:end,2:end-1)))./(dx(2:end,:)+dx(1:end-1,:)));
    yvalue=exp((dy(:,1:end-1).*log(phi.value(2:end-1,1:end-1))+...
        dy(:,2:end).*log(phi.value(2:end-1,2:end)))./(dy(:,2:end)+dy(:,1:end-1)));
    zvalue=[];
elseif (d == 3) || (d==3.2)
    Nx = phi.domain.dims(1);
    Ny = phi.domain.dims(2);
    Nz = phi.domain.dims(3);
    dx = repmat(phi.domain.cellsize.x, 1, Ny, Nz);
    dy = repmat(phi.domain.cellsize.y', Nx, 1, Nz);
    DZ = zeros(1,1,Nz+2);
    DZ(1,1,:) = phi.domain.cellsize.z;
    dz=repmat(DZ, Nx, Ny, 1);
    xvalue=exp((dx(1:end-1,:,:).*log(phi.value(1:end-1,2:end-1,2:end-1))+...
        dx(2:end,:,:).*log(phi.value(2:end,2:end-1,2:end-1)))./(dx(2:end,:,:)+dx(1:end-1,:,:)));
    yvalue=exp((dy(:,1:end-1,:).*log(phi.value(2:end-1,1:end-1,2:end-1))+...
        dy(:,2:end,:).*log(phi.value(2:end-1,2:end,2:end-1)))./(dy(:,1:end-1,:)+dy(:,2:end,:)));
    zvalue=exp((dz(:,:,1:end-1).*log(phi.value(2:end-1,2:end-1,1:end-1))+...
        dz(:,:,2:end).*log(phi.value(2:end-1,2:end-1,2:end)))./(dz(:,:,1:end-1)+dz(:,:,2:end)));
end
phiFaceAverage=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
