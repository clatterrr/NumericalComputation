function M=mass_matrixP1_2D(elements,areas,coeffs)
%FEM assembly of the mass matrix K representing the bilinear form 
%a(u,v) = (coeffs *u, v) for P1 functions u,v. 
%Here, coeffs represents a scalar diffusion coefficient specified as a collumn vector with number of entries equal to: 
%size(elements,1) - it represents an elementwise constant (P0) function.
%Note: the case of P1 coeffs requires a higher order quadrature (not implemented yet)
%If coeffs is not provided then coeffs=1 is assumed globally.


Xscalar=kron(ones(1,3),elements); Yscalar=kron(elements,ones(1,3)); 
if (nargin<3)
    Zmassmatrix=kron(areas,reshape((ones(3)+eye(3))/12,1,9)); 
else
    if numel(coeffs)==size(elements,1) %P0 coefficients
        Zmassmatrix=kron(areas.*coeffs,reshape((ones(3)+eye(3))/12,1,9)); 
    else %P1 coefficients
        M1=[6 2 2; 2 2 1; 2 1 2]/60;
        M2=M1([3,1,2],[3,1,2]);
        M3=M2([3,1,2],[3,1,2]);
            
        Zmassmatrix=kron(areas.*coeffs(elements(:,1)),reshape(M1,1,9)) ...
                   +kron(areas.*coeffs(elements(:,2)),reshape(M2,1,9)) ...
                   +kron(areas.*coeffs(elements(:,3)),reshape(M3,1,9));
    end
        
end

M=sparse(Xscalar,Yscalar,Zmassmatrix);