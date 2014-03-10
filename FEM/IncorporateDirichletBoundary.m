function [fnew Anew Mnew Unew dUnew] = IncorporateDirichletBoundary(A,M,U,dU,upperNodes,lowerNodes,upperPlate,lowerPlate)
lownode = sort(lowerNodes,'descend');
uppnode = sort(upperNodes,'descend');
allnode = sort([lowerNodes; upperNodes],'descend');
Anew = A;
Mnew = M;
Unew = U;
dUnew = dU;
%Removing rows:
for in=allnode
    i = 3*in; %z-coordinate.
    if i==size(Anew,1)
        Anew = Anew(1:(end-1),:);
        Mnew = Mnew(1:(end-1),:);
        Unew = Unew(1:(end-1));
        dUnew = dUnew(1:(end-1));
    else
        Anew = Anew([1:(i-1) (i+1):end],:);
        Mnew = Mnew([1:(i-1) (i+1):end],:);
        Unew = Unew([1:(i-1) (i+1):end]);
        dUnew = dUnew([1:(i-1) (i+1):end]);
    end
end

%Incorporating boundary conditions:
%Upper:
fnew = fnew -(Anew(:,upperNodes)+Mnew(:,upperNodes))*upperPlate;

%Lower:
fnew = fnew -(Anew(:,lowerNodes)+Mnew(:,lowerNodes))*lowerPlate;

%Removing columns:
for in=allnode
    i = 3*in;
    if i==size(Anew,2)
        Anew = Anew(:,1:(end-1));
        Mnew = Mnew(:,1:(end-1));
    else
        Anew = Anew(:,[1:(i-1) (i+1):end]);
        Mnew = Mnew(:,[1:(i-1) (i+1):end]);
    end
end
end