function [x_plates y_plates] = getSidePlates(p,tri)

boundary = unique(tri);
boundary_nodes = p(boundary,1:2);


x_plates = boundary(find(abs(boundary_nodes(:,1))==1));
y_plates = boundary(find(abs(boundary_nodes(:,2))==1));

end