function I = quadrature3D(p1, p2, p3, p4, Nq, g)
B = [(p2-p1) (p3-p1) (p4-p1)];
T = det(B)/6;

switch Nq
    case 1
        eta = 0.25*[1 1 1 1];
        weights = 1;
    case 4
        temp1 = 0.5854102;
        temp2 = 0.1381966;
        eta = [temp1 temp2 temp2 temp2; temp2 temp1 temp2 temp2; temp2 temp2 temp1 temp2; temp2 temp2 temp2 temp1];
        weights = 0.25*[1 1 1 1];
    case 5
        eta = 1/24*[6 6 6 6; 12 4 4 4; 4 12 4 4; 4 4 12 4; 4 4 4 12];
        weights = 1/20*[-16 9 9 9 9];
end

I = 0;
for i=1:Nq
    x= eta(i,1)*p1 + eta(i,2)*p2 + eta(i,3)*p3 + eta(i,4)*p4;
    I = I + T*weights(i)*feval(g,x);
end
end

