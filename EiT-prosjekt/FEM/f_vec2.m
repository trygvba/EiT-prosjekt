function [ f_vec ] = f_vec2( p, tri, force,epsilon,first_triangle )

nb=getNeumannBoundary2(tri,p,epsilon,first_triangle);



f_vec=zeros(3*length(p(:,1)),1);

N=size(nb,1);


if N<1
    f_vec=zeros(3*length(p(:,1)),1);
else
        for i=1:N
        c=nb(i,:);
        A(i)=3*planeIntegral3D(p(c(1),:),p(c(2),:),p(c(3),:));
        end

        A_tot=sum(A);

        Aa=A./A_tot;

        for i=1:N
            c=nb(i,:);
            f_vec(3*c)=f_vec(3*c)+Aa(i)*force*planeIntegral3D(p(c(1),:),p(c(2),:),p(c(3),:));
        end
end

end

