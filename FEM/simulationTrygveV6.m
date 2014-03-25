clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters:
%Polymer:
X = 15*10^(-6); %Length scale.
Ep = 0.8*10^9/X;;
vp = 0.46;
rhop = 950;

%Silver:
Es = 72.4*10^9/X;
vs = 0.37;
rhos = 10^4;

disp('Starting assembly')
%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

[p tri tetr] = loadGeo('spherewshell');
boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);

%--------------------------------------------
disp('Assembly done')
disp('Setting up for time integration')

%Parameters for time integration:
T0=0;                       %start time.
szU=size(A,1);              %dimension of our system.
steps=2000;                   %Number of time steps.
%U = zeros(szU,steps);
dt=10^(-9);                   %Temporal step size.
OLT=0.01;                   %Outer Layer Thickness.
impactzone=0.5;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=9*10^(6);                 %Frequency of upperplate.

%Getting out nodes on the Dirichlet boundary:
%Upper plate:
    bound_up = ballradius-0.025;
    [upperNodes, uz_up] = upperdirichletnodes(bound_up,p,zeros(szU,1),boundary);
    
%Lower plate:
    bound_low = -bound_up;
    [lowerNodes, uz_low] = lowerdirichletnodes(p,zeros(szU,1),bound_low,boundary);

%Setting initial displacement:
% U(3*upperNodes,1) = uz_up;
% U(3*lowerNodes,1) = uz_low;

%Setting initial velocity:

% v = zeros(szU,1);
% v(3*upperNodes) = -OLT*omega;
v=zeros(szU-length(upperNodes)-length(lowerNodes),1);

%------------------------------------------------

plateDisp = @(t) -OLT*sin(omega*t);
plateAcc = @(t) OLT*omega^2*sin(omega*t);

disp('Modifying matrices according to boundary conditions')
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

disp('Starting time integration')
%------------------------------------------------
%       TIME INTEGRATION (NEWMARK)
%------------------------------------------------
%Parameters for Paraview printing:
output_folder = 'paraview/animation';
title = 'testing';
NumberOfPics = 100;
n=1;
tic
for i=1:steps
    
    if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
        State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U);
        n = n+1;
    end
    t =T0+(i-1)*dt;
    F_cur = -Fmat_up*(uz_up+plateDisp(t))-Fmat_low*uz_low-plateAcc(t)*F_acc;
    Utemp_cur = K2*(Utemp_last+dt*v-0.25*dt^2*K1*Utemp_last+0.25*dt^2*(Mmod_inv*(F_cur+F_last)));
    utemp = putDirichletBack(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t));
    U = utemp;
    v = v+ 0.5*dt*(Mmod_inv*(F_cur+F_last))-0.5*dt*K1*(Utemp_cur+Utemp_last);
    Utemp_last = Utemp_cur;
    F_last = F_cur;
end
toc
disp('Done')