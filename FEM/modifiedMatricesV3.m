function [Anew Mnew Fmat_up Fmat_low F_acc] = modifiedMatricesV3(A,M,upperNodes,lowerNodes, x_plates, y_plates)
nodes_up = 3*upperNodes;
nodes_low = 3*lowerNodes;
nodes_x = 3*(x_plates-1)+1;
nodes_y = 3*(y_plates-1)+2;
allnode = sort([nodes_up; nodes_low; nodes_x; nodes_y],'descend');

Anew = A;
Mnew = M;
%Removing rows:
n=length(allnode);
for j=1:n
        i = allnode(j);
        Anew = Anew([1:(i-1) (i+1):end],:);
        Mnew = Mnew([1:(i-1) (i+1):end],:);
end

Fmat_up = Anew(:,3*upperNodes);
Fmat_low = Anew(:,3*lowerNodes);
F_acc = sum(Mnew(:,3*upperNodes),2);

%Removing columns:
for j=1:n
    i = allnode(j);
    Anew = Anew(:,[1:(i-1) (i+1):end]);
    Mnew = Mnew(:,[1:(i-1) (i+1):end]);

end
end
