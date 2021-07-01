%................................................................

% MATLAB codes for Finite Element Analysis
% problem9.m
% Bernoulli Beam in bending, under uniform load
% antonio ferreira 2008

% clear memory
clear all
clf
% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E=210000; A=100; I=2e8; EA=E*A; EI=E*I;
% generation of coordinates and connectivities
numberElements=6;
p1=3000*(1+cos(pi/4));
nodeCoordinates=[0 3000;
                 1500 3000;
                 3000 3000;
                 3000+1500*cos(pi/4) 1500;
                 p1 0;
                 p1+1500 0;
                 p1+3000 0];
xx=nodeCoordinates;
for i=1:numberElements; 
    elementNodes(i,1)=i; 
    elementNodes(i,2)=i+1;
end
numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
    % GDof: global number of degrees of freedom
GDof=3*numberNodes; 
U=zeros(GDof,1);
force=zeros(GDof,1);
stiffness=zeros(GDof,GDof); 
%force vector
force(10)=-10000;
force(12)=-10000;
force(17)=-5000000;
force(19)=5000000;

% calculation of the system stiffness matrix
% and force vector
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[ indice indice+numberNodes indice+2*numberNodes] 
  nn=length(indice);  
  xa=xx(indice(2))-xx(indice(1))
  ya=yy(indice(2))-yy(indice(1));  
  length_element=sqrt(xa*xa+ya*ya);
  cosa=xa/length_element;  
  sena=ya/length_element;
  ll=length_element;
   L= [cosa     0       sena    0       0       0;
       0        cosa    0       sena    0       0;
       -sena    0       cosa    0       0       0;
       0        -sena   0       cosa    0       0 ;
       0        0       0       0       1       0;
       0        0       0       0       0       1];

k1=[EA/ll   -EA/ll 0   0   0   0 ;
   -EA/ll  EA/ll   0   0   0   0 ;
0       0  12*EI/ll^3  -12*EI/ll^3  6*EI/ll^2   6*EI/ll^2;
0       0  -12*EI/ll^3 12*EI/ll^3  -6*EI/ll^2  -6*EI/ll^2;
0       0  6*EI/ll^2   -6*EI/ll^2  4*EI/ll     2*EI/ll;
0       0  6*EI/ll^2   -6*EI/ll^2  2*EI/ll     4*EI/ll]
    stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+L'*k1*L;
end 
% boundary conditions and solution
prescribedDof=[1 7 8 14 15 21];
activeDof=setdiff([1:GDof]', ...
    [prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;
% displacements
disp('Displacements')
jj=1:GDof; format
[jj' displacements]

U=displacements(1:3:GDof);   
U=displacements;   
clf
%desenho de malha original e deformadas
desenhoMalhas(nodeCoordinates+0.05*[U(1:numberNodes) U(numberNodes+1:2*numberNodes)],elementNodes,'L2','g.-');
desenhoMalhas(nodeCoordinates,elementNodes,'L2','b.-');


