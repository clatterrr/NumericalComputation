function [K,areas]=stiffness_matrixP1_2D(elements,coordinates,coeffs)
%FEM assembly of the stiffness matrix K representing the bilinear form 
%a(u,v) = (coeffs * grad u, grad v) for P1 functions u,v. 
%Here, coeffs represents a scalar diffusion coefficient specified as a collumn vector with number of entries equal to: 
%1) size(elements,1) - it represents an elementwise constant (P0) function,
%2) size(coordinates,1) - it represents an elementwise affine and globally continuous (P1) function. 
%If coeffs is not provided then coeffs=1 is assumed globally.

NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=3; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB
        coord(d,i,:)=coordinates(elements(:,i),d);
    end
end   
IP=[1/3;1/3];
[dphi,jac] = phider(coord,IP,'P1'); 
dphi = squeeze(dphi); 
areas=abs(squeeze(jac))/factorial(DIM);

if (nargin<3)
    Z=astam(areas',amtam(dphi,dphi));  
else
    if numel(coeffs)==size(coordinates,1)  %P1->P0 averaging
        coeffs=evaluate_average_point(elements,coeffs);
    end  
    Z=astam((areas.*coeffs)',amtam(dphi,dphi));
end
Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
X=permute(Y,[2 1 3]);
K=sparse(X(:),Y(:),Z(:));  