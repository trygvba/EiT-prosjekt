function [fnew Anew Mnew] = IncorporateDirichletBoundaryV2(A,M,upperNodes,lowerNodes,upperPlate,lowerPlate, accUpper)
lownode = sort(lowerNodes,'descend');
uppnode = sort(upperNodes,'descend');

allnode = sort([lowerNodes; upperNodes],'descend');

dim = length(A(1,:));
Anew = A;
Mnew = M;
%Removing rows:
n=length(allnode);
for j=1:n
    i = 3*allnode(j); %z-coordinate.
    %Removing columns from stiffness matrix:
    Anew(:,i) = sparse(dim,1);
    
    %Zeroing rows in A and M:
    Anew(i,:) = sparse(1,dim);
    Mnew(i,:) = sparse(1,dim);
    Mnew(i,i) = 1;  
end


%Incorporating boundary conditions:
fnew = sparse(size(Anew,1),1);

%Upper:
fnew = fnew -Anew(:,3*upperNodes)*upperPlate;

%Lower:
fnew = fnew -Anew(:,3*lowerNodes)*lowerPlate;

fnew = fnew -A(:,3*upperNodes)*upperPlate;

%Lower:
fnew = fnew -A(:,3*lowerNodes)*lowerPlate;


fnew(3*upperNodes) = accUpper*ones(length(upperNodes),1);
fnew(3*lowerNodes) = zeros(length(lowerNodes),1);
end

