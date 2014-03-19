function neu = getNeumannBoundary(tri,p,level)

neu=[];
for k=1:size(tri)
    pz1 = p(tri(k,1),3);
    pz2 = p(tri(k,2),3);
    pz3 = p(tri(k,3),3);
    if (pz1>=level)||(pz2>=level)||(pz3>=level)
        neu = [neu; tri(k,:)];
    end
end
end

