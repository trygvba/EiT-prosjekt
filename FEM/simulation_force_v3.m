clear all;
close all;
format long;
addpath(genpath('../Converters'));

%------------------------------------------
%       MAIN PARAMETERS:
%------------------------------------------
X = 15*10^(-6); %Length scale/physical radius.


%Declaration of parameters: 
%Polymer: 
Ep = 3*10^9*X;  
vp = 0.48; 
rhop = 950*X^3; 

%Silver: 
Es = 72.4*10^9*X; 
vs = 0.37; 
rhos = 10^4*X^3;
%--------------------------------------------

disp('Setting parameters done.')
disp('Loading geometry and assembling matrices.')

tic
%--------------------------------------------
%           ASSEMBLY:
%--------------------------------------------
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

[p tri tetr] = loadGeo('spherewshell');
boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);
%--------------------------------------------
toc
disp('Assembly done.')
disp('Setting up relevant parameters for time integration')
%--------------------------------------------
%Parameters for time integration:
T0=0;                       %Start time.
szU=size(A,1);              %Dimension of our system.
szP=szU/3;                  %Number of points.
steps=1000;                  %Number of time steps.
dt=10^(-8);                 %Temporal step size.
OLT=0.01;                   %Outer Layer Thickness.
impactzone=.05;             %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=100*pi;               %Frequency of upperplate.
howlow=-ballradius;         %Level of the lower plate.

f=8*10^(-2)/X;                %Force, needs to be calibrated.
epsilon=0.02;               %Epsilon layer of plate.


%Setting initial displacement and velocity:
v = zeros(szU,1);
U = zeros(szU,1);

%Upper plate position and acceleration as functions of time.
plateDisp = @(t) ballradius-OLT*sin(omega*t);
plateAcc = @(t) OLT*omega^2*sin(omega*t);

disp('Generating matrices to be used in time integration.')
K1 = M\A;
I = eye(szU);
K2 = (I+0.25*dt^2*K1)\I;
M_inv = M\I;
disp('Matrices are generated.')
%------------------------------------------------
%   TIME INTEGRATION (NEWMARK)
%------------------------------------------------
disp('Starting time integration.')
f_last = f_vec3(p,U,tri,-f,plateDisp(T0));
acc_last = M_inv*f_last;


output_folder = 'paraview/animation';
title = 'testing';
Number_of_pics = 100;
n=1;

 tic
 for i=1:(steps-1)    
     if i==1||(floor(Number_of_pics*i/steps)>floor(Number_of_pics*(i-1)/steps))
         State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U);
         n = n+1;
     end
     t = T0+i*dt;
     f_current = f_vec3(p,U,tri,-f,plateDisp(t));
     U = K2*(U+dt*v+0.25*dt^2*acc_last+0.25*dt^2*M_inv*f_current);
           

    [lower_nodes,uzi] = lowerdirichletnodes(p,U, howlow, boundary );
    U(3*lower_nodes)=uzi;
    acc_current = M_inv*(f_current-A*U);
    acc_current(3*lower_nodes) = 0;
    v = v+0.5*dt*(acc_last+acc_current);
    v(3*lower_nodes) = 0;
     acc_last = acc_current;
     f_last = f_current;     
 end
 toc
 disp('Done.')





