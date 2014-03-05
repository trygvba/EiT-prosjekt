function [ nodes ] = upperdirichletnodes( displacement, coords )

[N n]=size(coords);
nodes=[];
j=1;
for i=1:N
    
    if (sqrt(coords(i,1)^2 + coords(i,2)^2 +coords(i,3)^2)>=1.15) && (coords(i,3)>=displacement)
        nodes(j)=i;
        j=j+1;
    end
    
end

end
