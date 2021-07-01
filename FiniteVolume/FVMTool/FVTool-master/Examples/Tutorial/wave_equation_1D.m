% Written by Ali A. Eftekhari
% Last checked: June 2021
clc
c=@(x)(x.^2);
dc=@(x)(2*x);
Lx=1.0;
Nx=200;
dx=Lx/Nx;
m=createMesh1D(Nx, Lx);
x_face=m.facecenters.x;
x_cell=m.cellcenters.x;
BC=createBC(m);
BC.left.periodic=true;
BC.right.periodic=true;
u0=abs(sin(x_cell/Lx*10*pi));
u_old=createCellVariable(m, u0);
u_val=u_old;
dt=0.1;
c_face=createFaceVariable(m, 0.0);
c_face.xvalue=c(x_face);
dc_cell=createCellVariable(m, dc(x_cell));
Mconv=convectionUpwindTerm(c_face);
Ms=linearSourceTerm(dc_cell);
[Mbc, RHSbc]=boundaryCondition(BC);
for i=1:1000
  [Mt, RHSt]=transientTerm(u_old, dt, 1.0);
  M=Mt+Mconv-Ms+Mbc;
  RHS=RHSt+RHSbc;
  u_val=solvePDE(m, M, RHS);
  u_old=u_val;
  visualizeCells(u_val); drawnow;
end