clear rows level_size_P1 level_time_P1_CSV level_time_P1_RV

% add_paths;
addpath(genpath('3rd_party_elasticity_codes_for_testing'),'elasticity_3D');

levels=2; %maximum uniform refinement level

%homogeneous material parameters
young = 206900 ;                                     % Young's modulus E
poisson =  0.29 ;                                    % Poisson's ratio nu

lambda= young*poisson/((1+poisson)*(1-2*poisson)) ;  %Lamme first parameter
mu = young./(2*(1+poisson)) ;                        %Lamme second parameter

bulkC = young./(3*(1-2*poisson)) ;                    % bulk modulus K
shearC = mu ;                                         % shear modulus G (equal to Lamme second parameter)    
        
%coarse mesh
load mesh3D.mat; elements=elements3; clear elements3    

for level=0:levels    
    %uniform refinement
    if (level>0)
        [coordinates,elements]=refinement_uniform3D(coordinates,elements);
    end
    
    figure(4); visualize_mesh; 
    
    rows=numel(coordinates);
    level_size_P1(level+1)=rows; 
    
    %stiffness matrix assembly - method 1
    fprintf('technique of Cermak, Sysala and Valdman:\n')
    tic; 
    [Xi, WF] = quadrature_volume('P1');    % quadrature points and weights for volume integrals
    [HatP,DHatP1,DHatP2,DHatP3] = local_basis_volume('P1', Xi); % local basis functions and their derivatives for volume  integrals
    n_e=size(elements',2);     % number of elements
    n_q=length(WF);           % number of quadratic points
    n_int = n_e*n_q ;         % total number of integrations points
    shear =shearC*ones(1,n_int); bulk=bulkC*ones(1,n_int);
    [K,WEIGHT]=elastic_stiffness_matrix(elements',coordinates',shear,bulk,DHatP1,DHatP2,DHatP3,WF); 
    assembly_time=toc; 
    level_time_P1_CSV(level+1)=assembly_time; 
    
    rows=size(K,1);
    fprintf('level=%d, ', level);
    fprintf('time spent on K: %6.1e seconds, ',assembly_time(end));
    fprintf('rows of matrix =%d ',rows);
    fprintf('\n');  
    
    %stiffness matrix assembly - method 2
    fprintf('technique of Rahman and Valdman: \n')
    tic; 
    K2=stiffness_matrixP1_3D_elasticity(elements,coordinates,lambda,mu); 
    assembly_time=toc; 
    level_time_P1_RV(level+1)=assembly_time; 
    
    fprintf('level=%d, ', level);
    fprintf('time spent on K: %6.1e seconds, ',assembly_time(end));
    fprintf('rows of matrix =%d ',rows);
    fprintf('\n');  
    
    fprintf('-----------------------------------------------\n')
end

% output information
fprintf('\n')
for level=0:levels
    fprintf('%d ', level);
    fprintf('& ');
    fprintf('%d ', level_size_P1(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_P1_CSV(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_P1_RV(level+1));
    fprintf('\\\\');
    fprintf('\n');
end

rmpath(genpath('3rd_party_elasticity_codes_for_testing'),'elasticity_3D');
