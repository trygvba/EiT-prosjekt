function [A M] = MassAndStiffnessMatrix3D(tet,p,Cp,Cs,rhop,rhos)
N = size(tet,1);
dim = size(p,1);
Np = sum(tet(:,5)==62);

A =sparse(3*dim,3*dim);
M = sparse(3*dim,3*dim);

Mp = 1/120*[2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2];
I = eye(4);

for k=1:Np
    p1 = p(tet(k,1),:)';
    p2 = p(tet(k,2),:)';
    p3 = p(tet(k,3),:)';
    p4 = p(tet(k,4),:)';
    X = [(p2-p1) (p3-p1) (p4-p1)];
    T = det(X);
    V = abs(T/6);
    Coeff = [1 p1'; 1 p2'; 1 p3'; 1 p4'];
    
    for alpha=1:4
        ihat = tet(k,alpha);
        c_alpha = Coeff\I(:,alpha);
        i1 = 3*(ihat-1)+1;
        i2 = 3*(ihat-1)+2;
        i3 = 3*(ihat-1)+3;
        for beta = 1:4
            jhat = tet(k,beta);
            c_beta = Coeff\I(:,beta);
            j1 = 3*(jhat-1)+1;
            j2 = 3*(jhat-1)+2;
            j3 = 3*(jhat-1)+3;
            
            a1=[c_alpha(2) 0 0 c_alpha(3) c_alpha(4) 0];
            b1=[c_beta(2); 0;0; c_beta(3); c_beta(4); 0];
            a2=[0 c_alpha(3) 0 c_alpha(2) 0 c_alpha(4)];
            b2=[0;c_beta(3);0;c_beta(2);0;c_beta(4)];
            a3=[0 0 c_alpha(4) 0 c_alpha(2) c_alpha(3)];
            b3=[0;0;c_beta(4);0;c_beta(2);c_beta(3)];
            
            A(i1,j1) = A(i1,j1)+V*a1*Cp*b1;
            A(i1,j2) = A(i1,j2)+V*a1*Cp*b2;
            A(i1,j3) = A(i1,j3)+V*a1*Cp*b3;
            A(i2,j1) = A(i2,j1)+V*a2*Cp*b1;
            A(i2,j2) = A(i2,j2)+V*a2*Cp*b2;
            A(i2,j3) = A(i2,j3)+V*a2*Cp*b3;
            A(i3,j1) = A(i3,j1)+V*a3*Cp*b1;
            A(i3,j2) = A(i3,j2)+V*a3*Cp*b2;
            A(i3,j3) = A(i3,j3)+V*a3*Cp*b3;
            
            M(i1,j1) = M(i1,j1) + rhop*T*Mp(alpha,beta);
            M(i2,j2) = M(i2,j2) + rhop*T*Mp(alpha,beta);
            M(i3,j3) = M(i3,j3) + rhop*T*Mp(alpha,beta);
        end
    end
end


for k=(Np+1):N
    p1 = p(tet(k,1),:)';
    p2 = p(tet(k,2),:)';
    p3 = p(tet(k,3),:)';
    p4 = p(tet(k,4),:)';
    X = [(p2-p1) (p3-p1) (p4-p1)];
    T = det(X);
    V = abs(T/6);
    Coeff = [1 p1'; 1 p2'; 1 p3'; 1 p4'];
    
    for alpha=1:4
        ihat = tet(k,alpha);
        c_alpha = Coeff\I(:,alpha);
        i1 = 3*(ihat-1)+1;
        i2 = 3*(ihat-1)+2;
        i3 = 3*(ihat-1)+3;
        for beta = 1:4
            jhat = tet(k,beta);
            c_beta = Coeff\I(:,beta);
            j1 = 3*(jhat-1)+1;
            j2 = 3*(jhat-1)+2;
            j3 = 3*(jhat-1)+3;
            
            a1=[c_alpha(2) 0 0 c_alpha(3) c_alpha(4) 0];
            b1=[c_beta(2); 0;0; c_beta(3); c_beta(4); 0];
            a2=[0 c_alpha(3) 0 c_alpha(2) 0 c_alpha(4)];
            b2=[0;c_beta(3);0;c_beta(2);0;c_beta(4)];
            a3=[0 0 c_alpha(4) 0 c_alpha(2) c_alpha(3)];
            b3=[0;0;c_beta(4);0;c_beta(2);c_beta(3)];
            
            A(i1,j1) = A(i1,j1)+V*a1*Cs*b1;
            A(i1,j2) = A(i1,j2)+V*a1*Cs*b2;
            A(i1,j3) = A(i1,j3)+V*a1*Cs*b3;
            A(i2,j1) = A(i2,j1)+V*a2*Cs*b1;
            A(i2,j2) = A(i2,j2)+V*a2*Cs*b2;
            A(i2,j3) = A(i2,j3)+V*a2*Cs*b3;
            A(i3,j1) = A(i3,j1)+V*a3*Cs*b1;
            A(i3,j2) = A(i3,j2)+V*a3*Cs*b2;
            A(i3,j3) = A(i3,j3)+V*a3*Cs*b3;
            
            M(i1,j1) = M(i1,j1) + rhos*T*Mp(alpha,beta);
            M(i2,j2) = M(i2,j2) + rhos*T*Mp(alpha,beta);
            M(i3,j3) = M(i3,j3) + rhos*T*Mp(alpha,beta);
        end
    end
end

end

