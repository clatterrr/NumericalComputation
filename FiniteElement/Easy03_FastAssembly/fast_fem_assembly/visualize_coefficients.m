if (level<4)
    if numel(coeffs_stiffness_matrix)==size(elements,1)
        subplot(2,4,level+1); show_constant_scalar_frame(coeffs_stiffness_matrix,elements,coordinates); axis image; title(strcat('coeffs\_stiffness\_matrix, level=',num2str(level)));
    else
        subplot(2,4,level+1); show_nodal_scalar_frame(coeffs_stiffness_matrix,elements,coordinates); axis image; title(strcat('coeffs\_stiffness\_matrix, level=',num2str(level)));
    end
    
    if numel(coeffs_mass_matrix)==size(elements,1)        
        subplot(2,4,4+level+1); show_constant_scalar_frame(coeffs_mass_matrix,elements,coordinates); axis image; title(strcat('coeffs\_mass\_matrix, level=',num2str(level)));             
    else
        subplot(2,4,4+level+1); show_nodal_scalar_frame(coeffs_mass_matrix,elements,coordinates); axis image; title(strcat('coeffs\_mass\_matrix, level=',num2str(level)));             
    end
end