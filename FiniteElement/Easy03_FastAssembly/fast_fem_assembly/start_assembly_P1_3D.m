clear all

%addpath('.\library_vectorization')      %in windows
addpath('./library_vectorization/')      %in linux

levels=4; %maximum uniform refinement level for my comp with 8Gb memory is 7!

load mesh3D.mat; elements=elements3; clear elements3   %This replaces the mesh generation for older versions of Matlab, eg. R2008b

for level=0:levels    
    %uniform refinement
    if (level>0)
        [coordinates,elements]=refinement_uniform3D(coordinates,elements);
    end
    
    figure(2); visualize_mesh;
    
    %coeffs_stifness_matrix=sum(evaluate_average_point(elements,coordinates).^2,2); %P0 coefficient depending on a square of distance to the origin
    coeffs_stiffness_matrix=sum(coordinates.^2,2); %P1 coefficient depending on a square of distance to the origin
    w=[1 1/2 1/3]; coeffs_mass_matrix=coordinates*w'; %P1 coefficient in the form (x,y,z).w, where w is given vector
    %coeffs_mass_matrix=evaluate_average_point(elements,coeffs_mass_matrix); %recomputing P1 -> P0
    %coeffs_mass_matrix=ones(size(coeffs_mass_matrix));
    
    %stiffness and mass matrix assembly
    tic; [K,volumes]=stiffness_matrixP1_3D(elements,coordinates,coeffs_stiffness_matrix); time_stiffness_matrix(level+1)=toc; 
    tic; M=mass_matrixP1_3D(elements,volumes,coeffs_mass_matrix); time_mass_matrix(level+1)=toc;
    
    rows(level+1)=size(K,1);
    display_comparison
    
    if 0 %testing mass matrix 
        v1=coordinates(:,1)+coordinates(:,2)+coordinates(:,3); integral_value1(level+1)=v1'*M*v1;
        fprintf('   integral of (x+y/2+z/3)*(x+y+z)^2 dx dy dz over the unit cube = %f (the exact value is 11/4 = 2.75) \n', integral_value1(level+1) )
        
        v2=coordinates(:,1)+coordinates(:,2)-coordinates(:,3); integral_value2(level+1)=v2'*M*v2;   
        fprintf('   integral of (x+y/2+z/3)*(x+y-z)^2 dx dy dz over the unit cube = %f (the exact value is 5/9 = 0.55555555) \n', integral_value2(level+1) )     
    end
    
end