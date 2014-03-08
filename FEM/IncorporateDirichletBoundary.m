function [fnew Anew Mnew Unew] = IncorporateDirichletBoundary(A,M,U,upperNodes,lowerNodes,upperPlate,lowerPlate)
lownode = sort(lowerNodes,'descend');
uppnode = sort(upperNodes,'descend');
allnode = sort([lowerNodes upperNodes],'descend');
Anew = A;
Mnew = M;
Unew = U;
%Removing rows:
for in=allnode
    i = 3*in; %z-coordinate.
    if i==size(Anew,1)
        Anew = Anew(1:(end-1),:);
        Mnew = Mnew(1:(end-1),:);
        Unew = Unew(1:(end-1));
    else
        Anew = Anew([1:(i-1) (i+1):end],:);
        Mnew = Mnew([1:(i-1) (i+1):end],:);
        Unew = Unew([1:(i-1) (i+1):end]);
    end
end

%incorporating the upper plate:
dim = size(Anew,1);
fnew = sparse(dim,1);
for in=uppnode
    i = 3*in;
    fnew = fnew -upperPlate(in)*(Anew(:,i)+Mnew(:,i)); %Hvis det blir feil, vurder om Mnew ikke skal være med her!
end
%incorporating the lower plate:
for in=lownode
    i = 3*in;
    fnew = fnew -lowerPlate(in)*(Anew(:,i)+Mnew(:,i));
end

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