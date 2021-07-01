clear all

%addpath('.\library_vectorization')      %in windows
addpath('./library_vectorization/')      %in linux

levels=4; %maximum uniform refinement level
lambda=1; mu=1; %Lamme coeficients

load mesh3D.mat; elements=elements3; clear elements3    

for level=0:levels    
    %uniform refinement
    if (level>0)
        [coordinates,elements]=refinement_uniform3D(coordinates,elements);
    end
    
    figure(4); visualize_mesh;
        
    %stiffness matrix assembly
    tic; [K,areas]=stiffness_matrixP1_3D_elasticity(elements,coordinates,lambda,mu); time1(level+1)=toc; 
    rows(level+1)=size(K,1);
       
    fprintf('level=%d, ', level);
    fprintf('time spent on stifness matrix K: %f seconds, ',time1(level+1));
    fprintf('size of square matrice K =%d ',rows(level+1));
    fprintf('\n');
end


