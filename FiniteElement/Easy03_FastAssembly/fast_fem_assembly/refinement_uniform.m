function [coordinates,elements3,dirichlet]=refinement_uniform(coordinates,elements3,dirichlet) 
%function: [coordinates,elements3,dirichlet]=refinement_uniform(coordinates,elements3,dirichlet) 
%requires: getEdges, symetrizeMatrix  
%uniform refinement of a 2D triangulation 

%uniform refinement   
[element2edges, edge2nodes]=getEdges(elements3);    
nodes2edge=sparse(edge2nodes(:,1),edge2nodes(:,2),1:size(edge2nodes,1),size(coordinates,1),size(coordinates,1));
nodes2edge=symetrizeMatrix(nodes2edge); 
    
%elements on uniformly refined mesh
elements3_internal=element2edges+size(coordinates,1);
elements3_refin1= [elements3(:,1) elements3_internal(:,3) elements3_internal(:,2)];
elements3_refin2= [elements3(:,2) elements3_internal(:,1) elements3_internal(:,3)];
elements3_refin3= [elements3(:,3) elements3_internal(:,2) elements3_internal(:,1)];    
elements3=[elements3_internal; elements3_refin1; elements3_refin2; elements3_refin3];  

if (nargin==3)
%dirichlet edges of uniformly refined mesh
dirichlet_edges=diag(nodes2edge(dirichlet(:,1),dirichlet(:,2)));
dirichlet=[dirichlet(:,1) dirichlet_edges+size(coordinates,1); dirichlet_edges+size(coordinates,1) dirichlet(:,2)];
end


%coordinates of uniformly refined mesh
coordinates_internal=(coordinates(edge2nodes(:,1),:)+coordinates(edge2nodes(:,2),:))/2;
coordinates=[coordinates; coordinates_internal];   

function A_sym = symetrizeMatrix(A)
[i,j,k]=find(A);
W=sparse([i; j], [j; i], ones(size(k,1)*2,1));
A_help=sparse([i; j], [j; i], [k; k]);
[i,j,k]=find(A_help);
[i,j,kk]=find(W);
A_sym=sparse(i,j,(kk.^(-1)).*k); %Now Kantennr_sym is a symetric form of Kantennr