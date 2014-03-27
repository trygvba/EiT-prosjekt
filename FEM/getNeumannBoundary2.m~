function neu = getNeumannBoundary2(tri,p,epsilon,first_triangle)

neu=[];
ft=mean(p(first_triangle,3));

for k=1:size(tri)
    pz1 = p(tri(k,1),3);
    pz2 = p(tri(k,2),3);
    pz3 = p(tri(k,3),3);
    if ((abs(pz1- ft) <=epsilon))||((abs(pz2- ft) <=epsilon))||((abs(pz1- ft) <=epsilon))
        neu = [neu; tri(k,:)];
    end
    
end
end

