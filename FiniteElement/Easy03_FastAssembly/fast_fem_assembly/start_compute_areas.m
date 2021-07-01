%maximal level of mesh refinement for my computer with 8Gb memory is 11 
%(it creates a triangular mesh with 25 milion elements)! 
levels=9; 

coordinates=[0 0; 1 0; 2 0; 2 1; 1 1; 0 1;0 2; 1 2];
elements=[1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7];

fprintf('mesh refinement - level: ');
for level=0:levels
    if (level>0)
        [coordinates,elements]=refinement_uniform(coordinates,elements);   
        fprintf('%d ', level);
    end

    figure(6); visualize_mesh;
end
fprintf('\n');
fprintf('the finest mesh ready: ');
fprintf('number of elements=%d, ',size(elements,1));
fprintf('number of nodes=%d \n',size(coordinates,1));
fprintf('\n');

if 1
    fprintf('vectorized version, ');
    tic
    a=coordinates(elements(:,1),1);
    b=coordinates(elements(:,2),1);
    c=coordinates(elements(:,3),1);
    d=coordinates(elements(:,1),2);
    e=coordinates(elements(:,2),2);
    f=coordinates(elements(:,3),2);
    areas=abs(b.*f - c.*e - a.*f + c.*d + a.*e - b.*d)/2; 
    time_vectorized=toc; 
else
    fprintf('vectorized version 2, ');
    tic
    a=coordinates(elements(:,1),1);
    b=coordinates(elements(:,2),1);
    c=coordinates(elements(:,3),1);
    d=coordinates(elements(:,1),2);
    e=coordinates(elements(:,2),2);
    f=coordinates(elements(:,3),2);
    ca=c-a;
    cb=c-b;
    fd=f-d;
    fe=f-e;
    areas2=abs(ca.*fe - cb.*fd)/2; 
    time_vectorized=toc;
end
fprintf('time spent: %3.2f seconds \n', time_vectorized); 

fprintf('serial version, ');
tic
areas=zeros(size(elements,1),1);
for level=1:size(elements,1)
    a=coordinates(elements(level,1),1);
    b=coordinates(elements(level,2),1);
    c=coordinates(elements(level,3),1);
    d=coordinates(elements(level,1),2);
    e=coordinates(elements(level,2),2);
    f=coordinates(elements(level,3),2);
    areas(level)=abs(a*e + b*f + c*d - a*f - b*d - c*e)/2;
end    
time_serial=toc;
fprintf('time spent: %3.2f seconds \n', time_serial);
fprintf('\n');

fprintf('the vectorized code is %3.1f times faster than the serial code! \n', time_serial/time_vectorized);





