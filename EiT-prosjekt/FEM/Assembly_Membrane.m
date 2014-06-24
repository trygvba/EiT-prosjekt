function [p, tri, tetr, szU, K1, K2, Amod, Mmod_inv, vel,Fmat_up, Fmat_low, F_acc, lowerNodes, uz_low,upperNodes, uz_up,omega, dt, x_plates, y_plates] = Assembly_Membrane(GeoName,Phys_groups, E,v,rho,N, granularity)

disp('Starting assembly')
%-----------ASSEMBLY:------------------------
C=zeros(6,6*length(Phys_groups));

for i = 1:length(Phys_groups)
    C(:,6*i-5:6*i) = StressMatrix(E(i),v(i));
end

[p, tri, tetr] = loadGeo(GeoName);
boundary = unique(tri);
[A, M] = MassAndStiffnessMatrixGeneral(tetr,p,C, rho, Phys_groups);%,Cs,rhop,rhos);
spy(M)
%--------------------------------------------
disp('Scaling time and frequency')
szU=size(A,1);              %dimension of our system.
maxz=max(p(:,3));     %Total radius of the ball with outer shell.
lambdap = E(1)*(v(1))/((1+max(v))*(1-2*max(v)));
mup = E(1)/(2*(1+max(v)));
harmonicK = pi/(4*maxz); %Standing wave number with node at x such that kx = pi/2
omega= N*harmonicK*sqrt((lambdap + 2*mup)/rho(1)); % Frequency of upperplate by w = kv. Coefficient should be an odd k-multiple for damping.
dt = granularity/omega; % Timestep granularity
%Getting out nodes on the Dirichlet boundary:

disp('Modifying matrices according to boundary conditions')

%Upper plate:
    bound_up = maxz;
    [upperNodes, uz_up] = upperdirichletnodes(bound_up,p,zeros(szU,1),boundary);
    
%Lower plate:
    bound_low = -1*bound_up;
    [lowerNodes, uz_low] = lowerdirichletnodes(p,zeros(szU,1),bound_low,boundary);
  
%Side plates:
    [x_plates, y_plates] = getSidePlates(p,tri);


    

vel=zeros(szU-length(upperNodes)-length(lowerNodes)-length(x_plates)-length(y_plates),1);


%Getting out modified matrices:
[Amod Mmod Fmat_up Fmat_low F_acc] = modifiedMatricesV3(A,M,upperNodes,lowerNodes, x_plates, y_plates);
I = eye(size(Amod,1));
Mmod_inv = Mmod\I;
K1 = Mmod_inv*Amod;
K2 = (I+(0.25*dt^2)*K1)\I;

disp('Assembly done')
end