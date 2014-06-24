function volume = tetvol(p,tet,k)
    p1 = p(tet(k,1),:)';
    p2 = p(tet(k,2),:)';
    p3 = p(tet(k,3),:)';
    p4 = p(tet(k,4),:)';
    X = [(p2-p1) (p3-p1) (p4-p1)];
    T = det(X);
    volume = abs(T/6);
end