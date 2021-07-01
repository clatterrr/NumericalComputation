% =========================================================================
%
%  This program triggers an assembly test for a 3D elastoplastic body. It 
%  is considered perfectly plastic model with the Drucker-Prager yield 
%  criterion and the strip-footing benchmark. The main aim is to compare
%  the assembling time for the elastic and tangent stiffness matrices. 
%  The tangent stiffness matrices are computed in each time step. 
%  One can set optionally 4 types of finite elements,
%  levels of mesh density and many other parameters.
%
% ======================================================================
%

%
% The main input data 
%

  elem_type='Q2'; % the type of finite elements; available choices:
                  % 'P1', 'P2', 'Q1', 'Q2'
  level=0;        % a nonnegative integer defining mesh density
      
  % values of elastic material parameters
  young = 1e7;                           % Young's modulus
  poisson =  0.48;                       % Poisson's ratio
  shear = young/(2*(1+poisson)) ;        % shear modulus
  bulk = young/(3*(1-2*poisson)) ;       % bulk modulus    
  
  % values of plastic material parematers
  c0 = 450 ;                                 % cohesion
  phi = pi/9;                                % frictional angle  
  eta = 6*sin(phi)/(sqrt(3)*(3+sin(phi)));   % inner approximation
  c = 6*c0*cos(phi)/(sqrt(3)*(3+sin(phi)));  % inner approximation
       
%
% Mesh generation
%

  % geometrical parameters (choose only integers)
  size_xy = 10;      % size of the body in direction x and y 
  size_z = 1;        % size of the body in z-direction  
                      
  % the mesh generation depending prescribed finite elements        
  switch(elem_type)
    case 'P1'
        [COORD,ELEM,SURF,DIRICHLET,Q]=mesh_P1(level,size_xy,size_z);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD,ELEM,SURF,DIRICHLET,Q]=mesh_P2(level,size_xy,size_z);
        fprintf('P2 elements: \n')
    case 'Q1'
        [COORD,ELEM,SURF,DIRICHLET,Q]=mesh_Q1(level,size_xy,size_z);
        fprintf('Q1 elements: \n')
    case 'Q2'
        [COORD,ELEM,SURF,DIRICHLET,Q]=mesh_Q2(level,size_xy,size_z);
        fprintf('Q2 elements: \n')
    otherwise
          disp('bad choice of element type');
  end  
  Q_ND=DIRICHLET(2,:)>0;

%
% Data from the reference element
%
  
  % quadrature points and weights for volume integration
  [Xi, WF] = quadrature_volume(elem_type);
  
  % local basis functions and their derivatives for volume
  [HatP,DHatP1,DHatP2,DHatP3] = local_basis_volume(elem_type, Xi);

%
% Number of nodes, elements and integration points + print
%
  n_n=size(COORD,2);          % number of nodes
  n_unknown=length(COORD(Q)); % number of unknowns
  n_e=size(ELEM,2);           % number of elements
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 
  % 
  fprintf('number of nodes =%d ',n_n);  
  fprintf('\n');   
  fprintf('number of unknowns =%d ',n_unknown);
  fprintf('\n');   
  fprintf('number of elements =%d ',n_e);
  fprintf('\n');   
  fprintf('number of integration points =%d ',n_int);
  fprintf('\n');   

%
% Assembling of the elastic stiffness matrix
%
  
  % values of elastic material parameters at integration points       
  shear=shear*ones(1,n_int);
  bulk=bulk*ones(1,n_int); 
   
  % stiffness matrix assembly and the assembly time
  tic;     
  [K_elast,B,WEIGHT,iD,jD,D_elast]=elastic_stiffness_matrix(ELEM,COORD,...
                            shear,bulk,DHatP1,DHatP2,DHatP3,WF);  
  assembly_elast_time=toc; 
  fprintf('step =%d ',1);
  fprintf('\n');   
  fprintf('  time spent on K_elast:  %6.1e seconds, ',assembly_elast_time);
  fprintf('\n');   

%
% Plastic material parameters
%

  % values of plastic material parematers at integration points
  eta=eta*ones(1,n_int);
  c=c*ones(1,n_int); 
   
%
% Loading process
%    

  % initial load increment and factor
  d_zeta=1/1000;   % load increment
  d_zeta_min=d_zeta/1300;
  d_zeta_old=d_zeta;
  zeta=0;          % load factor
  zeta_old=zeta;
  zeta_max=1;      % maximal value of the load factor  
  
  % elastic initial displacements
  Ud = -d_zeta*DIRICHLET;
  f = - K_elast*Ud(:);
  U_it = Ud;
  U_it(Q) = K_elast(Q,Q)\f(Q) ;  
  
  % other initialization
  dU = zeros(3,n_n) ;      % Newton's increment (in displacement)
  U = zeros(3,n_n) ; 
  U_old = -U_it ; 
  F = zeros(3,n_n) ;       % vector of internal forces
  E = zeros(6,n_int);      % strain tensors at integration points
  Ep_old = zeros(6,n_int); % plastic strain tensors at integration points
  zeta_hist=zeros(1,1000);
  pressure_old=0;
  pressure_hist=zeros(1,1000);
  
  % storage of assembly time in dependence on plastic integration points
  assembly=zeros(20000,2);
  assembly_step=0;   
  AUX=reshape(1:6*n_int,6,n_int);
  
  % loop over load steps 
  step=1;
  fprintf('step =%d ',step);
  while true  
    
    % updated load factor  
    zeta=zeta_old+d_zeta;  
    fprintf('load factor =%d ',zeta);  
    fprintf('load increment =%d ',d_zeta);  
    fprintf('pressure =%d ',pressure_old);  
    fprintf('\n');             
     
    % Newton's solver (the semismooth Newton method)
    n_it=25;           % maximal number of Newton iterations
    it=0;              % iteration number
    while it<n_it       
      
      it=it+1;   
        
      % consistent tangent stiffness matrix
      tic
      E(:) = B*U_it(:) ;   % strain at integration points
      [S,DS,IND_p]=constitutive_problem(E,Ep_old,shear,bulk,eta,c);
      vD = repmat(WEIGHT,36,1).*DS ; 
      D_p = sparse( iD(:),jD(:),vD(:), 6*n_int,6*n_int ) ;   
      K_tangent = K_elast+B'*(D_p-D_elast)*B;     
      assembly_time=toc;
      
      % measuring assembly dependance on plastic integration points
      n_plast=length(WEIGHT(IND_p));
      assembly_step=assembly_step+1;
      assembly(assembly_step,:)=[n_plast assembly_time];   
      fprintf('  time spent on K_tangent: %6.1e seconds, ',assembly_time);
       
      % vector of internal forces
      F(:) = B'*reshape(repmat(WEIGHT,6,1).*S, 6*n_int,1) ;
      
      % Newton's increment
      dU(Q) = K_tangent(Q,Q)\(-F(Q)); 
                   
      % next iteration
      U_new= U_it + dU ;

      % stopping criterion 
      q1 = sqrt( dU(:)'*K_elast*dU(:) ) ;
      q2 = sqrt(  U_it(:)'*K_elast*U_it(:)  ) ;
      q3 = sqrt( U_new(:)'*K_elast*U_new(:) ) ;
      criterion = q1/(q2+q3);
      if isnan(criterion)
          it=n_it;
          break
      end
      fprintf('  stopping criterion=%6.1e  ',criterion); 
      fprintf('\n');  
      
      % update of unknown arrays
      U_it=U_new;                                                     
            
      % test on the stopping criterion
      if  criterion < 1e-12
          break
      end    
          
    end%  Newton
     
    if criterion<1e-10  % successful convergence of the Newton method
      U_old=U;
      U=U_it;
      E(:) = B*U(:) ;
      [S,DS,IND_p,lambda,Ep]=constitutive_problem(E,Ep_old,shear,bulk,eta,c); 
      Ep_old=Ep;
      zeta_old=zeta;
      d_zeta_old=d_zeta;
      step=step+1; 
      zeta_hist(step)=zeta;
      % normalized mean pressure on the footing
      pressure_array=transformation(S(2,:),ELEM,WEIGHT);
      pressure=-mean(full(pressure_array(Q_ND)))/c0;
      pressure_hist(step)=pressure;
      if (pressure-pressure_old<0.05)&&(criterion<1e-12)
          d_zeta=2*d_zeta;
      end
      pressure_old=pressure;
    else
      warning('The Newton solver does not converge.')
      d_zeta=d_zeta/2; % decrease of the load increment
    end
    
    % initialization for the next iteration
    U_it=d_zeta*(U-U_old)/d_zeta_old+U;   
    
    % stopping criteria for the loading process
    
    if zeta_old>=zeta_max
        warning('Maximal load factor was achieved.')
        break
    end
    
    if d_zeta<d_zeta_min
          warning('Too small load increments.')
          break
    end    
    
    fprintf('step =%d ',step);
 
  end % loading process
  

%  
% visualization of selected results
%
      
  % mesh visualization
  if level<2
    draw_mesh(COORD,SURF,elem_type)     
  end
  
  % load path
  figure
  semilogx(zeta_hist(1:step),pressure_hist(1:step));
  xlabel('settlement'); ylabel('normalized pressure');
     
  % total displacements + deformed shape
  U_total = sqrt(U(1,:).^2 + U(2,:).^2 + U(3,:).^2);
  draw_quantity(COORD,SURF,U,U_total,elem_type,size_xy,size_z)
  colorbar off; colorbar('location','south')  

%   fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
%   print('-painters','-dpdf',strcat('displacement_',num2str(step)))    
    
  % plastic multipliers
  Ep_node = transformation(sqrt(sum(Ep.*Ep)),ELEM,WEIGHT); 
  draw_quantity(COORD,SURF,0*U,zeros(size(Ep_node))+Ep_node,elem_type,size_xy,size_z)
 
  % dependance of the assembly time on plastic integration points
  figure
  x=assembly(1:assembly_step,1); y=assembly(1:assembly_step,2); 
  x_ext=n_int; y_ext_meas=assembly_elast_time;
  x_r=x(x>0); y_r=y(x>0);  
 
  X = [ones(length(x_r),1) x_r];
  b = X\y_r;  % linear regression, y ~ b(0)+b(1)*x
  x_e=[x; x_ext];      
  y_e=[ones(length(x_e),1) x_e]*b;
  y_ext=[1 x_ext]*b;
  plot(x,y,'bx',x_e,y_e,'r-',x_ext,y_ext,'ro',x_ext,y_ext_meas,'bo');
  legend('plasticity - assembly times','plasticity - linear fit',...
  'plasticity - extrapolation','elasticity - assembly time','location','northwest')    
      


