function [fnew Anew Mnew Unew dUnew] = IncorporateDirichletBoundaryV3(A,M,U,upperNodes,lowerNodes,upperPlate,lowerPlate, plateAcc)
lownode = sort(lowerNodes,'descend');
uppnode = sort(upperNodes,'descend');

allnode = sort([lowerNodes; upperNodes],'descend');

Anew = A;
Mnew = M;
Unew = U;
%Removing rows:
n=length(allnode);
for j=1:n
    i = 3*allnode(j); %z-coordinate.
    Anew = Anew([1:(i-1) (i+1):end],:);
    Mnew = Mnew([1:(i-1) (i+1):end],:);
    Unew = Unew([1:(i-1) (i+1):end]);
end

%Incorporating boundary conditions:

fnew = sparse(size(Anew,1),1);

%Upper:
fnew = fnew -Anew(:,3*upperNodes)*upperPlate-Mnew(:,3*upperNodes)*ones(length(upperNodes),1)*plateAcc;

%Lower:
fnew = fnew -Anew(:,3*lowerNodes)*lowerPlate;

%Removing columns:
for j=1:n
    i = 3*allnode(j);
    Anew = Anew(:,[1:(i-1) (i+1):end]);
    Mnew = Mnew(:,[1:(i-1) (i+1):end]);

end
end