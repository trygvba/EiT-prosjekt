function [ bdcoords ] = tricoords( tri, p )


bdcoords=[p(tri(:,1),1:3) p(tri(:,2),1:3) p(tri(:,3),1:3)] ;

end

