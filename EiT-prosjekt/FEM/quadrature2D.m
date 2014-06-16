function I = quadrature2D(p1, p2, p3, Nq, g)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B =[(p2-p1) (p3-p1)];
T = 0.5*det(B);
    
switch Nq
    case 1
        eta = 1/3*[1 1 1];
        weights = 1;
    case 3
        eta = 0.5*[1 1 0; 1 0 1; 0 1 1];
        weights = 1/3*[1 1 1];
    case 4
        eta = 1/15*[5 5 5; 9 3 3; 3 9 3; 3 3 9];
        weights = 1/48*[-27 25 25 25];
end

I=0;
for i=1:Nq
    x = p1*eta(i,1)+p2*eta(i,2)+p3*eta(i,3);
    I = I + T*weights(i)*feval(g,x);
end
end

