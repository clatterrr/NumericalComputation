if (level<3)
    subplot(1,3,level+1); show_mesh(elements,coordinates); axis image; title(strcat('mesh, level=',num2str(level)));
end