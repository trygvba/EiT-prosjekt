function [p, tri, tetr, szU, K1, K2, Amod, Mmod_inv, v,Fmat_up, Fmat_low, F_acc, lowerNodes, uz_low,upperNodes, uz_up,omega, dt] = HomogenousAssembly(GeoName,Ep,vp,rhop,N)

disp('Starting assembly')
%-----------ASSEMBLY:------------------------
Cp = StressMatrix(Ep,vp);

[p tri tetr] = loadGeo(GeoName);
boundary = unique(tri);
[A M] = HomogenousMaterial(tetr,p,Cp, rhop);%,Cs,rhop,rhos);

%--------------------------------------------
disp('Scaling time and frequency')
szU=size(A,1);              %dimension of our system.
maxz=max(p(:,3));     %Total radius of the ball with outer shell.
lambdap = Ep*vp/((1+vp)*(1-2*vp));
mup = Ep/(2*(1+vp));
harmonicK = pi/(4*maxz); %Standing wave number with node at x such that kx = pi/2
omega= N*harmonicK*sqrt((lambdap + 2*mup)/rhop); % Frequency of upperplate by w = kv. Coefficient should be an odd k-multiple for damping.
dt = 0.05/omega; % Timestep granularity
%Getting out nodes on the Dirichlet boundary:

disp('Modifying matrices according to boundary conditions')

%Upper plate:
    bound_up = maxz;
    [upperNodes, uz_up] = upperdirichletnodes(bound_up,p,zeros(szU,1),boundary);
    
%Lower plate:
    bound_low = -1*bound_up;
    [lowerNodes, uz_low] = lowerdirichletnodes(p,zeros(szU,1),bound_low,boundary);
    

    

v=zeros(szU-length(upperNodes)-length(lowerNodes),1);


%Getting out modified matrices:
[Amod Mmod Fmat_up Fmat_low F_acc] = modifiedMatricesV2(A,M,upperNodes,lowerNodes);
K1 = Mmod\Amod;
I = eye(size(Amod,1));
K2 = (I+0.25*dt^2*K1)\I;
Mmod_inv = Mmod\I;

%Utemp_last = v;
F_last = -Fmat_up*uz_up-Fmat_low*uz_low;
Utemp_last = Amod\F_last;
u0 = putDirichletBack(Utemp_last,lowerNodes,upperNodes,uz_low,uz_up);
U = u0;
disp('Assembly done')
end