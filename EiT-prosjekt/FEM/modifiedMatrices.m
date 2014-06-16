function [Anew Mnew] = modifiedMatrices(A,M,upperNodes,lowerNodes)
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
    %zeroing columns from stiffness matrix:
    Anew(:,i) = 0;
    
    %Zeroing rows in A and M:
    Anew(i,:) = 0;
    Mnew(i,:) = 0;
    Mnew(i,i) = 1;  
end
end


