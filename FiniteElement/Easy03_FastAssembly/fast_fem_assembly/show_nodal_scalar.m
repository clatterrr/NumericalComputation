function show_nodal_scalar(nodalValue,elements,coordinates,nodalDisplacement)
switch nargin, 
    case 3,
        nodalDisplacement=0*coordinates;
    case {0, 1, 2}
        fprintf('missing parameters')
end


X=coordinates(:,1)+nodalDisplacement(:,1);
Y=coordinates(:,2)+nodalDisplacement(:,2);

fill3(X(elements)',Y(elements)',nodalValue(elements)',nodalValue(elements)','FaceColor','interp','LineStyle','none');
set(gcf,'renderer','zbuffer');