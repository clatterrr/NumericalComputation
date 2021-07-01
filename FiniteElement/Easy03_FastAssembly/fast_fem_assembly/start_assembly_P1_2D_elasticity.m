clear all

%addpath('.\library_vectorization')     %in windows
addpath('./library_vectorization/')     %in linux

levels=8; %maximum uniform refinement level

lambda=1; mu=1; %Lamme coeficients

create_2D_mesh; %mesh for testing

for level=0:levels    
    %uniform refinement
    if (level>0)
        [coordinates,elements,dirichlet]=refinement_uniform(coordinates,elements,dirichlet);
    end
    
    figure(3); visualize_mesh;
    
    %stiffness and mass matrix assembly
    tic; [K areas]=stiffness_matrixP1_2D_elasticity(elements,coordinates,lambda,mu); time_stiffness_matrix(level+1)=toc; 
    tic; M=mass_matrixP1_2D_elasticity(elements,areas); time_mass_matrix(level+1)=toc;
    
    rows(level+1)=size(K,1); 
    display_comparison
end


