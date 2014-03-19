clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters:
%Polymer:
Ep = 10; %0.8*10^9;
vp = 0.46;
rhop = 950;

%Silver:
Es = 72.4/0.8*Ep; %72.4*10^9;
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
steps=20;                   %Number of time steps.
U = zeros(szU,steps);
dt=10^(-10);                   %Temporal step size.
OLT=0.05;                   %Outer Layer Thickness.
impactzone=0.5;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=1;                 %Frequency of upperplate.




%Setting initial velocity:
v = zeros(szU,1);
v(3*upperNodes) = -OLT*omega;
%------------------------------------------------

plateDisp = @(t) -OLT*sin(omega*t);
plateAcc = @(t) OLT*omega^2*sin(omega*t);

plateForce=@(pD) = %double derivative of platedisp?







%------------------------------------------------
%   TIME INTEGRATION (NEWMARK)
%------------------------------------------------

for i=2:steps
    t = T0+i*dt;
    F_cur = -Fmat_up*(uz_up+plateDisp(t))-Fmat_low*uz_low;
    F_cur(3*upperNodes) = plateAcc(t);
    F_cur(3*lowerNodes) = 0;
    
    U(:,i) = K2*(U(:,i-1)+dt*v-0.25*dt^2*K1*U(:,i-1)+0.25*dt^2*Mmod\(F_cur+F_last));
    v = v+0.5*dt*Mmod\(F_cur+F_last)-0.5*dt*K1*(U(:,i-1)+U(:,i));
    
    F_last = F_cur;
end
     


%------------POST PROCESSING:-----------------
output_folder = 'paraview/animation';
title = 'testing';

for n=1:steps
    State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,n));
end


