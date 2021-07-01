%................................................................

% MATLAB codes for Finite Element Analysis
% problem21.m
% free vibrations of laminated plates
% See reference:
% K. M. Liew, Journal of Sound and Vibration, 
% Solving the vibration of thick symmetric laminates 
% by Reissner/Mindlin plate theory and the p-Ritz method, Vol. 198,
% Number 3, Pages 343-360, 1996

% antonio ferreira 2008

% clear memory
clear all;colordef white;clf

% materials
thickness=0.001;h=thickness;kapa=pi*pi/12;
rho=1;I=thickness^3/12;
% symbolic computation 
syms phi pi

% liew material
e2=1;e1=40*e2;g23=0.5*e2;g13=0.6*e2;g12=g13;
miu12=0.25;miu21=miu12*e2/e1;factor=1-miu12*miu21;

% angles for laminate
alfas=[0,pi/2,0];% 3 layers
% upper and lower coordinates
z(1)=-(h/2);z(2)=-(h/2)+h/3;z(3)=-z(2);z(4)=-z(1);

% [Q] in 0� orientation 
qbarra(1,1,1)=e1/factor;
qbarra(1,2,1)=miu21*e1/factor;
qbarra(2,1,1)=miu12*e2/factor;
qbarra(2,2,1)=e2/factor;
qbarra(3,3,1)=g12;
qbarra(4,4,1)=kapa*g23;
qbarra(5,5,1)=kapa*g13;

% transformation matrix
T=[cos(phi)^2,sin(phi)^2,-sin(2*phi),0,0;...
   sin(phi)^2,cos(phi)^2,sin(2*phi),0,0;...
   sin(phi)*cos(phi),-sin(phi)*cos(phi),cos(phi)^2-sin(phi)^2,0,0;...
   0,0,0,cos(phi),sin(phi);...
   0,0,0,-sin(phi),cos(phi)];

% [Q] in structural axes
qBarra=T*qbarra*T.';

for s=1:size(alfas,2)
    for i=1:5
        for j=1:5
           QQbarra(i,j,s)=subs(qBarra(i,j,1),phi,alfas(s));
       end
   end
   Qbarra=double(QQbarra);
end
Q=Qbarra; 

%______________________________________________
Astiff(5,5)=0;Bstiff(5,5)=0;Fstiff(5,5)=0;Istiff(5,5)=0;
for k=1:size(alfas,2)
        for i=1:3
        for j=1:3
        Astiff(i,j)=Astiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
        Bstiff(i,j)=Bstiff(i,j)+Q(i,j,k)*(z(k+1)^2-z(k)^2)/2;
        Fstiff(i,j)=Fstiff(i,j)+Q(i,j,k)*(z(k+1)^3-z(k)^3)/3;
        end
        end

            for i=4:5
            for j=4:5
            Istiff(i,j)=Istiff(i,j)+Q(i,j,k)*(z(k+1)-z(k));
            end
            end
end
pi=double(pi); % come back to numeric computation

% constitutive matrices
CMembranaMembrana=Astiff(1:3,1:3);
CMembranaFlexao0=Bstiff(1:3,1:3);
CFlexao0Flexao0=Fstiff(1:3,1:3);
CCorte0Corte0=Istiff(4:5,4:5);

% load
P = -1;

%Mesh generation
L  = 1;    
numberElementsX=10;
numberElementsY=10;
numberElements=numberElementsX*numberElementsY;

[nodeCoordinates, elementNodes] = ...
    rectangularMesh(2*L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
axis off
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=5*numberNodes; 

% stiffness and mass matrices      
stiffness=formStiffnessMatrixMindlinQ45laminated5dof...
    (GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,CMembranaMembrana,...
CMembranaFlexao0,CFlexao0Flexao0,CCorte0Corte0);

[mass]=...
    formMassMatrixMindlinQ4laminated5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,rho,thickness,I);

% boundary conditions 
[prescribedDof,activeDof,fixedNodeW]=...
    EssentialBC5dof('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);

% eigenproblem: free vibrations
numberOfModes=12;

[V,D] = eig(stiffness(activeDof,activeDof),...
    mass(activeDof,activeDof)); 
% Liew, p-Ritz
    D0=e2*h^3/12/(1-miu12*miu21);
    D = diag(sqrt(D)*L*L/pi/pi*sqrt(rho*h/D0));
    [D,ii] = sort(D); ii = ii(1:numberOfModes); 
    VV = V(:,ii);
    activeDofW=setdiff([1:numberNodes]',[fixedNodeW]);
    NNN=size(activeDofW);
    
    VVV(1:numberNodes,1:12)=0;
    for i=1:numberOfModes
        VVV(activeDofW,i)=VV(1:NNN,i);
    end
 
NN=numberNodes;N=sqrt(NN);
x=linspace(-L,L,numberElementsX+1);
y=linspace(-L,L,numberElementsY+1);


% drawing Eigenmodes
drawEigenmodes2D(x,y,VVV,NN,N,D)
