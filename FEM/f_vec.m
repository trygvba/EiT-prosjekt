function [ f_vec ] = f_vec( p, tri, force,displacement )

nb=getNeumannBoundary(tri,p,displacement);

f_vec=zeros(3*length(p(:,1)),1);

for i=1:length(nb(:,1))
    c=nb(i,:);
    f_vec(3*c)=f_vec(3*c)+force*planeIntegral3D(p(c(1),:),p(c(2),:),p(c(3),:));
end


end

