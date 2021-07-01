clear all

%addpath('.\library_vectorization')      %in windows
addpath('./library_vectorization/')      %in linux

levels=7; 
lambda=1; mu=1; %Lamme coeficients

load hexmesh.mat   %8 meshes (only partially nested) up to level 7, provided by TIEN DAT NGO (EPFL Lausanne), created in LS-PrePost

for level=0:levels   
    %uniform refinement
    elements=hex{level+1}.elements;
    coordinates=hex{level+1}.coordinates;

    figure(5); visualize_mesh;
     
    %stiffness matrix assembly
    tic; [K,areas]=stiffness_matrixQ1_3D_elasticity(elements,coordinates,lambda,mu); time1(level+1)=toc; 

    rows(level+1)=size(K,1);
       
    fprintf('level=%d, ', level);
    fprintf('time spent on K %6.1e seconds, ',time1(level+1));
    %fprintf('time spent on mass matrix M=%%10.3e, ',time2(level+1));
    fprintf('size of square matrix K=%d ',rows(level+1));
    fprintf('\n');
end

%creating hexahedral mesh on a unit cube - not working yet!
% NE_X=2;
% NE_Y=1;
% NE_Z=1;
% [X Y Z] = meshgrid(0:(1/NE_X):1,0:(1/NE_Y):1,0:(1/NE_Z):1);
% NN=(NE_X+1)*(NE_Y+1)*(NE_Z+1);
% coordinates=[reshape(X,NN,1,1) reshape(Y,NN,1,1) reshape(Z,NN,1,1)];
% element=[1 2 3+NE_X 2+NE_X [1 2 3+NE_X 2+NE_X]+(NE_X+1)*(NE_Y+1)];
% elements=element;
% for i=1:NE_X-1
%     elements=[elements; element+i];
% end


