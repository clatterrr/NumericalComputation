
function [coord,elem,surf,dirichlet,Q]=mesh_P1(level,size_xy)

% =========================================================================
%
%  This function creates triangular mesh for P1 elements
%
%  input data:
%    level   - an integer defining a density of a uniform mesh
%    size_xy - size of the body in directions x and y (integer)
%    body=(0,size_xy)  x (0,size_xy)
%
%  output data:
%    coord     - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%                number of nodes
%    elem      - array containing numbers of nodes defining each element, 
%                size(elem)=(3,n_e), n_e = number of elements
%    surf      - array containing numbers of nodes defining each surface element, 
%                size(surf)=(2,n_s), n_s = number of surface elements
%    dirichlet - array indicating the nodes where the nonhomogenous
%                Dirichlet boundary condition is considered,
%                size(dirichlet)=(2,n_n)
%    Q         - logical array indicating the nodes where the homogeneous
%                Dirichlet boundary condition is considered, size(Q)=(2,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_xy*2^level;      % number of segments in x direction
  N_y = N_x;                  % number of segments in y direction
 
  % 
  n_e = 2*N_x*N_y    ;        % total number of elements

%
% C - 2D auxilliary array that contains node numbers and that is important 
% for the mesh construction.
%
  C=reshape(1:(N_x+1)*(N_y+1),N_x+1,N_y+1);
 
%
% coordinates of nodes
%
  % coordinates in directions x and y
  coord_x=linspace(0,size_xy,N_x+1);
  coord_y=linspace(0,size_xy,N_y+1);
  
  % long 1D arrays containing coordinates of all nodes in x,y directions
  c_x=repmat(coord_x,1,N_y+1);     
  c_y=repmat(kron(coord_y,ones(1,N_x+1)),1);  
         
  % the required array of coordinates, size(coord)=(2,n_n)    
  coord=[c_x; c_y] ;  

% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0], V2 -> [1 0], V3 -> [1 1], V4 -> [0 1]
  %  V1,...,V4 are logical 2D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(N_x+1,N_y+1);
  V1(1:N_x,1:N_y)=1;
  %
  V2=false(N_x+1,N_y+1);
  V2(2:(N_x+1),1:N_y)=1;
  %
  V3=false(N_x+1,N_y+1);
  V3(2:(N_x+1),2:(N_y+1))=1;
  %
  V4=false(N_x+1,N_y+1);
  V4(1:N_x,2:(N_y+1))=1; 

  % used division of a rectangle into 2 triangles:   
  %   V1 V2 V4
  %   V2 V3 V4
  % size(aux_elem)=(2*3,n_e/2)
  aux_elem=[C(V1)'; C(V2)'; C(V4)';
            C(V2)'; C(V3)'; C(V4)' ];
        
  % the array elem, size(elem)=(3,n_e)
  elem=reshape(aux_elem,3,n_e);     
 
%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C
  % 
  % and logical 2D arrays V1_s, V2_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit line:
  %   V1_s -> [0 0], V2_s -> [1 0]      
  
  % Edge 1: y=0 (the bottom of the body)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,1);  
  V1_s=false(N_x+1,1);  V1_s(1:N_x    ,1)=1;
  V2_s=false(N_x+1,1);  V2_s(2:(N_x+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf1=reshape(aux_surf,2,N_x);       
  
  % Edge 2: x=size_xy (the right hand side of the body)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(end,:);  
  V1_s=false(N_y+1,1);  V1_s(1:N_y    ,1)=1;
  V2_s=false(N_y+1,1);  V2_s(2:(N_y+1),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf2=reshape(aux_surf,2,N_y);       
  
  % Edge 3: y=size_xy (the top of the body)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,end);  
  V1_s=false(N_x+1,1);  V1_s(1:N_x    ,1)=1;
  V2_s=false(N_x+1,1);  V2_s(2:(N_x+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf3=reshape(aux_surf,2,N_x);       
  
  % Edge 4: x=0 (the left hand side of the body)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(1,:);  
  V1_s=false(N_y+1,1);  V1_s(1:N_y    ,1)=1;
  V2_s=false(N_y+1,1);  V2_s(2:(N_y+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf4=reshape(aux_surf,2,N_y);     
   
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4] ;
  
%
% Boundary conditions
%   

  % array indicating the nodes with non-homogen. Dirichlet boundary cond.
  dirichlet = zeros(size(coord));
  dirichlet(2,(coord(2,:)==size_xy)&(coord(1,:)<=1.0001)) = 1;     
  
  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  Q(2,(coord(2,:)==size_xy)&(coord(1,:)<=1.0001)) = 0;      
  Q(1,(coord(1,:)==size_xy)) = 0;   
  
end
