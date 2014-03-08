function [nodes,uzi] = lowerdirichletnodes( p,u, howlow, boundary )

temp = [boundary ((p(boundary,3)+u(3*boundary))<=howlow)];
nodes = temp(find(temp(:,2)),1);
uzi = howlow*ones(length(nodes),1)-p(nodes,3);


end

