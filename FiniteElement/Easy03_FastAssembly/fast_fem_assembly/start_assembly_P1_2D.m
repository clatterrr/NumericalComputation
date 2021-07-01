%clear all

%addpath('.\library_vectorization')     %in windows
addpath('./library_vectorization/')     %in linux

levels=9; %maximal level of refinement, e.g. levels=11 for L-shape geometry creates FEM matrices with 12.5 million of rows 

create_2D_mesh; %mesh for testing

for level=0:levels    
    %uniform refinement
    if (level>0)
        [coordinates,elements,dirichlet]=refinement_uniform(coordinates,elements,dirichlet);
        %save(strcat(strcat('meshTriangularLevel',num2str(level)),'.txt'),'coordinates','elements','-ASCII','-tabs','-double');
    end
    
    %coeffs_stifness_matrix=sum(evaluate_average_point(elements,coordinates).^2,2); %P0 coefficient depending on a square of distance to the origin
    coeffs_stiffness_matrix=sum(coordinates.^2,2); %P1 coefficient depending on a square of distance to the origin
    g=[1 1/2]; coeffs_mass_matrix=coordinates*g'; %P1 coefficient in the form (x,y,z).w, where w is given vector
    %coeffs_mass_matrix=evaluate_average_point(elements,coeffs_mass_matrix); %recomputing P1 -> P0
    %coeffs_mass_matrix=ones(size(coeffs_mass_matrix));
    
    figure(1); visualize_coefficients;
    figure(11); visualize_mesh;
    
    %stiffness and mass matrix assembly
    tic; [K,areas]=stiffness_matrixP1_2D(elements,coordinates, coeffs_stiffness_matrix); time_stiffness_matrix(level+1)=toc; 
    tic; M=mass_matrixP1_2D(elements,areas,coeffs_mass_matrix); time_mass_matrix(level+1)=toc;
    
    
    if 0  %testing mass matrix 
        v=coordinates(:,1)+coordinates(:,2);    
        integral_value(level+1)=v'*M*v;
        fprintf('   integral of (x+y/2)*(x+y)^2 dx dy over the unit square = %f (the exact value is 1.125) \n', integral_value(level+1) )
    end
    
    
    rows(level+1)=size(K,1); 
    display_comparison
end




