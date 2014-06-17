function [ f_vec ] = f_vec3( p,U, tri, force,plateDisp)
p_new = pupdate(p,U);       %Setting absolute coordinates.
nb=getNeumannBoundary(tri,p_new,plateDisp); %Getting out which triangles are on the contact surface.



f_vec=zeros(3*length(p(:,1)),1);    %Initializing the loading vector.

N=size(nb,1);                       %Number of triangles.


if N<1
    f_vec=zeros(3*length(p(:,1)),1);
else
        %Iterate over each triangle
        for i=1:N
        c=nb(i,:);      %Get out relevant node indices.
        A(i)=3*planeIntegral3D(p_new(c(1),:),p_new(c(2),:),p_new(c(3),:));  %Area of each triangle on the surface.
        end

        A_tot=sum(A);   %Total area of contact surface.

        Aa=A./A_tot;    %Creating area weights.

        for i=1:N
            c=nb(i,:);  %Get out relevant node indices.
            %Adding contribution by (area weight)*(downward force)*(nodal
            %basis function integrated over triangle)
            f_vec(3*c)=f_vec(3*c)+Aa(i)*force*A(i)/3;
        end
end

end



