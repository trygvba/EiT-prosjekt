function I = planeIntegral3D(p1,p2,p3)

B = cross(p2-p1,p3-p1);
I = 1/6*sqrt(B*B');
end

