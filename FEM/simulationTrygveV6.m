clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters:
%Polymer:
Ep = 0.8*10^9;
vp = 0.46;
rhop = 950;

%Silver:
Es = 72.4*10^9;
vs = 0.37;
rhos = 10^4;


%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

[p tri tetr] = loadGeo('spherewshell');
boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);

%--------------------------------------------
%Parameters for time integration:
T0=0;                       %start time.
szU=size(A,1);              %dimension of our system.
steps=100;                   %Number of time steps.
U = zeros(szU,steps);
dt=10^(-6);                   %Temporal step size.
OLT=0.1;                   %Outer Layer Thickness.
impactzone=0.5;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=(pi/2)/(steps*dt);                 %Frequency of upperplate.

%Getting out nodes on the Dirichlet boundary:
%Upper plate:
    bound_up = ballradius-0.025;
    [upperNodes, uz_up] = upperdirichletnodes(bound_up,p,zeros(szU,1),boundary);
    
%Lower plate:
    bound_low = -bound_up;
    [lowerNodes, uz_low] = lowerdirichletnodes(p,zeros(szU,1),bound_low,boundary);

%Setting initial displacement:
U(3*upperNodes,1) = uz_up;
U(3*lowerNodes,1) = uz_low;

%Setting initial velocity:

% v = zeros(szU,1);
% v(3*upperNodes) = -OLT*omega;
v=zeros(szU-length(upperNodes)-length(lowerNodes),1);

%------------------------------------------------

plateDisp = @(t) -OLT*sin(omega*t);
plateAcc = @(t) OLT*omega^2*sin(omega*t);

%Getting out modified matrices:
[Amod Mmod Fmat_up Fmat_low F_acc] = modifiedMatricesV2(A,M,upperNodes,lowerNodes);
K1 = Mmod\Amod;
I = eye(size(Amod,1));
K2 = (I+0.25*dt^2*K1)\I;

%Utemp_last = v;
F_last = -Fmat_up*uz_up-Fmat_low*uz_low;
Utemp_last = Amod\F_last;
u0 = putDirichletBack(Utemp_last,lowerNodes,upperNodes,uz_low,uz_up);
U(:,1) = u0;
%------------------------------------------------
%       TIME INTEGRATION (NEWMARK)
%------------------------------------------------
for i=2:steps
    t =T0+(i-1)*dt;
    F_cur = -Fmat_up*(uz_up+plateDisp(t))-Fmat_low*uz_low-plateAcc(t)*F_acc;
    Utemp_cur = K2*(Utemp_last+dt*v-0.25*dt^2*K1*Utemp_last+0.25*dt^2*(Mmod\(F_cur+F_last)));
    utemp = putDirichletBack(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t));
    U(:,i) = utemp;
    v = v+ 0.5*dt*(Mmod\(F_cur+F_last))-0.5*dt*K1*(Utemp_cur+Utemp_last);
    Utemp_last = Utemp_cur;
    F_last = F_cur;
end





output_folder = 'paraview/animation';
title = 'testing';
for n=1:steps
    State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,n));
end
