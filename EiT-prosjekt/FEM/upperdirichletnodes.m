function [ nodes, uzi ] = upperdirichletnodes( displacement, p, u, boundary )

temp = [boundary ((p(boundary,3)+u(3*boundary))>=displacement)];
nodes = temp(find(temp(:,2)),1);
uzi = displacement*ones(length(nodes),1)-p(nodes,3);

end
