function show_constant_scalar_frame(constantValue,elements,coordinates,nodalDisplacement)
    
if (nargin==4)     
    if norm(nodalDisplacement)==0
       factor=1;
    else
       factor =10^(-round(log10(max(max(nodalDisplacement)))));
    end   
else
    nodalDisplacement=zeros(size(coordinates));
end

X=coordinates(:,1)+nodalDisplacement(:,1);
Y=coordinates(:,2)+nodalDisplacement(:,2);
fill3(X(elements)',Y(elements)',kron(ones(1,size(elements,2)),constantValue)',kron(ones(1,size(elements,2)),constantValue)','FaceColor','interp','Linewidth',0.1);



end
