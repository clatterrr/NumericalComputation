function [element2faces, face2edges, face2nodes]=getFaces(elements)
%function: [element2faces, face2nodes]=getFaces(elements)
%requires: deleterepeatedrows
%generates faces of (tetrahedral 3D) triangulation defined in elements
%elements is matrix, whose rows contain numbers of its element nodes 
%face2edges returns faces numbers of each tetrahedral element
%face2nodes returns three node numbers of each face
%example: [element2faces, face2nodes]=...
%getFaces([1 3 4 5; 7 4 3 5; 5 7 6 4; 6 8 7 4; 2 1 4 5; 2 4 6 5])

% %extracts sets of faces as sets of their nodes (vertices)
% nodes=zeros(size(elements,1)*4,3);
% nodes(1:4:end,:)=elements(:,[2 3 4]);
% nodes(2:4:end,:)=elements(:,[3 4 1]);
% nodes(3:4:end,:)=elements(:,[4 1 2]);
% nodes(4:4:end,:)=elements(:,[1 2 3]);
% 
% %repeated sets of nodes (joint faces) are eliminated 
% [face2nodes,element2faces]=deleterepeatedrows(nodes);
% element2faces=reshape(element2faces,size(elements,2),size(elements,1))';

%extracts sets of faces edges
[element2edges, edge2nodes]=getEdges(elements);
edges=zeros(size(elements,1)*4,3);
edges(1:4:end,:)=element2edges(:,[4 5 6]);
edges(2:4:end,:)=element2edges(:,[2 3 6]);
edges(3:4:end,:)=element2edges(:,[1 3 5]);
edges(4:4:end,:)=element2edges(:,[1 2 4]);

%repeated sets of nodes (joint faces) are eliminated 
[face2edges,element2faces]=deleterepeatedrows(edges);
element2faces=reshape(element2faces,size(elements,2),size(elements,1))';

face2nodes=[edge2nodes(face2edges(:,1),:), edge2nodes(face2edges(:,2),:) edge2nodes(face2edges(:,3),:)];




