rows(level+1)=size(K,1);   
fprintf('level=%d, ', level);
fprintf('time spent on K: %6.1e seconds, ',time_stiffness_matrix(level+1));
fprintf('time spent on M: %6.1e seconds, ',time_mass_matrix(level+1));
fprintf('size of square matrices K,M =%d ',rows(level+1));
fprintf('\n');