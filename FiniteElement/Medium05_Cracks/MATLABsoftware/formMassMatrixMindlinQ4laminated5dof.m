%................................................................

function [M]=...
    formMassMatrixMindlinQ4laminated5dof(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,rho,thickness,I)

% computation of mass matrix
% for Mindlin plate element

M=zeros(GDof);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
 
% cycle for element
for e=1:numberElements       
  % indice : nodal condofectivities for each element
  indice=elementNodes(e,:);           
  ndof=length(indice);
  
  % cycle for Gauss point
  for q=1:size(gaussWeights,1)                      
    GaussPoint=gaussLocations(q,:);                                                     
    xi=GaussPoint(1);
    eta=GaussPoint(2);

% shape functions and derivatives
    [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);

% Jacobian matrix, inverse of Jacobian, 
% derivatives w.r.t. x,y    
    [Jacob,invJacobian,XYderivatives]=...
        Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
    
% [B] matrix bending
    B_b=zeros(3,5*ndof);
    B_b(1,ndof+1:2*ndof)        = XYderivatives(:,1)';  
    B_b(2,2*ndof+1:3*ndof)      = XYderivatives(:,2)';
    B_b(3,ndof+1:2*ndof)        = XYderivatives(:,2)';  
    B_b(3,2*ndof+1:3*ndof)      = XYderivatives(:,1)';
% [B] matrix membrane
    B_m=zeros(3,5*ndof);
    B_m(1,3*ndof+1:4*ndof)      = XYderivatives(:,1)';  
    B_m(2,4*ndof+1:5*ndof)      = XYderivatives(:,2)';
    B_m(3,3*ndof+1:4*ndof)      = XYderivatives(:,2)';  
    B_m(3,4*ndof+1:5*ndof)      = XYderivatives(:,1)';
    
% mass matrix 
    M(indice,indice)=M(indice,indice)+...
        shapeFunction*shapeFunction'*thickness*rho*...
        gaussWeights(q)*det(Jacob);
    
    M(indice+numberNodes,indice+numberNodes)=...
    M(indice+numberNodes,indice+numberNodes)+...
        shapeFunction*shapeFunction'*I*rho*...
        gaussWeights(q)*det(Jacob);
    
    M(indice+2*numberNodes,indice+2*numberNodes)=...
        M(indice+2*numberNodes,indice+2*numberNodes)+...
        shapeFunction*shapeFunction'*I*rho*...
        gaussWeights(q)*det(Jacob);
    
    M(indice+3*numberNodes,indice+3*numberNodes)=...
        M(indice+3*numberNodes,indice+3*numberNodes)+...
        shapeFunction*shapeFunction'*thickness*rho*...
        gaussWeights(q)*det(Jacob);
    
    M(indice+4*numberNodes,indice+4*numberNodes)=...
        M(indice+4*numberNodes,indice+4*numberNodes)+...
        shapeFunction*shapeFunction'*thickness*rho*...
        gaussWeights(q)*det(Jacob);
    end  % Gauss point
end    % element

