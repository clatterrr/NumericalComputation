function  elements2midpoint=evaluate_average_point(elements,coordinates,which_elements_number)
        
dummy=coordinates(elements(:,1),:);
for i=2:size(elements,2)
    dummy=dummy+coordinates(elements(:,i),:);
end
elements2midpoint=dummy/size(elements,2);

if nargin==3
   elements2midpoint=elements2midpoint(which_elements_number,:);
end

   
end
