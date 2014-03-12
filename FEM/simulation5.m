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


%Parameters for time integration:
T0=0;                       %start time.
szU=size(A,1);              %dimension of our system.
steps=20;                   %Number of time steps.
U = zeros(szU,steps);
dt=0.001;                   %Temporal step size.
OLT=0.05;                   %Outer Layer Thickness.
impactzone=0.5;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=2*pi;                 %Frequency of upperplate.

%Getting out nodes on the Dirichlet boundary:
%Upper plate:
    bound_up = ballradius-impactzone*OLT;
    [upperNodes, uz_up] = upperdirichletnodes(bound_up,p,zeros(szU,1),boundary);
    
%Lower plate:
    bound_low = -bound_up;
    [lowerNodes, uz_low] = lowerdirichletnodes(p,zeros(szU,1),bound_low,boundary);

U(3*upperNodes,1) = uz_up;
U(3*lowerNodes,1) = uz_low;

output_folder = 'paraview/animation';
title = 'testing';
State_to_vtk(output_folder,title,1,szU,tetr(:,1:4),p,U(:,1));