function [Anew Mnew Fmat_up Fmat_low F_acc] = modifiedMatricesV2(A,M,upperNodes,lowerNodes)
allnode = sort([lowerNodes; upperNodes],'descend');

Anew = A;
Mnew = M;
%Removing rows:
n=length(allnode);
for j=1:n
    i = 3*allnode(j); %z-coordinate.
    Anew = Anew([1:(i-1) (i+1):end],:);
    Mnew = Mnew([1:(i-1) (i+1):end],:);
end

Fmat_up = Anew(:,3*upperNodes);
Fmat_low = Anew(:,3*lowerNodes);
F_acc = sum(Mnew(:,3*upperNodes),2);

%Removing columns:
for j=1:n
    i = 3*allnode(j);
    Anew = Anew(:,[1:(i-1) (i+1):end]);
    Mnew = Mnew(:,[1:(i-1) (i+1):end]);

end
end

