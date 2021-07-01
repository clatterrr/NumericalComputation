function draw_quantity(coord,elem,U,Q_node,elem_type,size_xy)

% =========================================================================
%
%  This function depicts prescribed nodal quantity
%
%  input data:
%    coord     - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%                number of nodes
%    elem      - array containing numbers of nodes defining each element, 
%                size(elem)=(n_p,n_e), n_e = number of elements
%    U         - nodal displacements, size(U)=(2,n_n) to catch deformed shape
%                if the deformed shape is not required then set 0*U
%    Q_node    - prescribed nodal quantity, size(Q_node)=(1,n_n)
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%    size_xy   - size of the body in x and y direction (integer)
%    body=(0,size_xy)x(0,size_xy)
%
% ======================================================================
%

  figure;
  hold on;
  
  % visualization of the quantity
  if strcmp(elem_type,'P1')||strcmp(elem_type,'P2')
    s = patch('Faces',elem(1:3,:)','Vertices',coord'+U',...
        'FaceVertexCData',Q_node','FaceColor','interp','EdgeColor','none'); 
  else
    s = patch('Faces',elem(1:4,:)','Vertices',coord'+U',...
        'FaceVertexCData',Q_node','FaceColor','interp','EdgeColor','none'); 
  end
  if ~isOctave  
    alpha(s,.5);
  end
  colorbar;
  
 % undeformed shape of the body
  plot([0,size_xy],[0,0])
  plot([0,size_xy],[size_xy,size_xy])
  plot([0,0],[0,size_xy])
  plot([size_xy,size_xy],[0,size_xy]) 
  
  %
  box on
  view(2);      % standard view in 2D
  axis equal;   % real ratios
  hold off;
  axis off;
end