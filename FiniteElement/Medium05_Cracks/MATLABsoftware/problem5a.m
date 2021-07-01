%................................................................

% MATLAB codes for Finite Element Analysis
% problem4.m
% ref: D. Logan, A first couse in the finite element method,
% third Edition, half-structure
% antonio ferreira 2008

% clear memory
clear all
% E; modulus of elasticity
% A: area of cross section
% L: length of bar
E=70000;   A=300;    EA=E*A;
% generation of coordinates and connectivities
elementNodes=[ 1 2;1 3;2 3;2 4;1 4;3 4];
nodeCoordinates=[ 0 0;0 3000;3000 0;3000 3000];
numberElements=size(elementNodes,1);
numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
% for structure:
    % displacements: displacement vector
    % force : force vector
    % stiffness: stiffness matrix
U=zeros(2*numberNodes,1);
force=zeros(2*numberNodes,1);
stiffness=zeros(2*numberNodes,2*numberNodes); 
% applied load at node 2
force(4)=-50000;
force(8)=-50000;
% calculation of the system stiffness matrix
for e=1:numberElements; 
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));
  length_element=sqrt(xa*xa+ya*ya);
  C=xa/length_element;
  S=ya/length_element;   
    k1=EA/length_element*...
        [C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;
        -C*C -C*S C*C C*S;-C*S -S*S C*S S*S];    
  stiffness(elementDof,elementDof)=...
      stiffness(elementDof,elementDof)+k1;
end 
% boundary conditions and solution
prescribedDof=[1 2 5 7]';
activeDof=setdiff([1:2*numberNodes]', ...
    [prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(2*numberNodes,1);
displacements(activeDof)=U;
us=1:2:2*numberNodes-1;
vs=2:2:2*numberNodes;
% displacements
disp('Displacements')
jj=1:2*numberNodes; format
[jj' displacements]

% drawing displacements

format long
figure
L=xx(2)-xx(1);
%L=node(2,1)-node(1,1);
XX=displacements(us);YY=displacements(vs);
dispNorm=max(sqrt(XX.^2+YY.^2));
scaleFact=2*dispNorm;
clf
hold on
drawingMesh(nodeCoordinates+scaleFact*[XX YY],elementNodes,'L2','k.-');
drawingMesh(nodeCoordinates,elementNodes,'L2','k.--');

% reactions
F=stiffness*displacements;
reactions=F(prescribedDof);
disp('reactions');format
[prescribedDof reactions]

% stresses at elements
for e=1:numberElements                          
  indice=elementNodes(e,:);
  elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));
  length_element=sqrt(xa*xa+ya*ya);
  C=xa/length_element;
  S=ya/length_element;   
  sigma(e)=E/length_element* ...
      [-C  -S C S]*displacements(elementDof); 
end    
disp('stresses')
sigma'