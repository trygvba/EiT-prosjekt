function [ nodes ] = lowerdirichletnodes( coords, howlow )

N=length(coords(:,1));
nodes=[];
j=1;
for i=1:N
    
    if (sqrt(coords(i,1)^2 + coords(i,2)^2 +coords(i,3)^2)>=1.18) && (coords(i,3)<= howlow)
        nodes(j)=i;
        j=j+1;
    end
    
end



end

