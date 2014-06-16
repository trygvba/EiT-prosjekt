%-----------------------------------------------------|
%   Script for running a FEM simulation on a sphere   |
%   with a thin coating.                              |
%-----------------------------------------------------|

clear all;
close all;
format long;

%---------------------------
%Declaration of parameters:
%---------------------------

X = 15*10^(-6)/2;     %Length scale.
%Polymer: 
Ep = 2*10^9*X;        %Young's modulus.
vp = 0.3;             %Poisson's ration.
rhop = 1.02*10^3*X^3; %Mass density. 

%Silver: 
Es = 80*10^9*X; 
vs = 0.35; 
rhos = 10.5*10^3*X^3;

%-------------------
%   ASSEMBLY:
%-------------------

%Creating the stress-strain matrices:
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

%Loading the geometry created in gmsh.
[p tri tetr] = loadGeo('spherewshell');
%Getting out which nodes lie on the boundary.
boundary = unique(tri);

%Constructing mass- and stiffness matrices:
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);

%-----------------------------------------
%       TIME INTEGRATION:
%-----------------------------------------

%Setting up parameters for time integration:
szU = size(A,1); %number of degrees of freedom.
epsilon =0.05;  %Layer thickness for Neumann approach.
ballradius = max(p(:,3));   %Finding radius of total sphere (with shell).
    
dt=10^(-9);     %Temporal step size.
steps = 1000;   %Number of time steps.

U = zeros(szu,steps);   %Initialising output data structure.
v = U(:,1);             %Initialising velocity vector.

%Setting up iteration matrices:
K1 = M\A;
K2 = (eye(szU) + beta*dt^2*K1)\eye(szU);

%Initialising force vector:
f_last = force(0,epsilon,p,U);      %Function that calculates force applied
                                    %at time t and distributes it along
                                    %relevant boundary found using epsilon,
                                    %p and U(:,i).

%Start of actual time integration:
for i=1:(steps-1)
    t = i*dt;
    f_current = force(t,epsilon,p,U(:,i)); %Force at current step
    
    %Updating displacement:
    U(:,i+1) = K2*(U(:,i) + dt*v - (1-2*beta)/2*dt^2*K1*U(:,i) ...
             + 0.5*dt^2*(M\((1-2*beta)*f_last+2*beta*f_current)));
    %Updating velocity:
    v = v + 0.5*dt*K1*(U(:,i)+U(:,i+1)) ...
      + 0.5*dt*(M\((1-2*beta)*f_last+2*beta*f_current));
    
    %Finding nodes in contact with the lower plate:
    [nodes, uz_low]=lowerdirichletnodes(p,U(:,i+1), -ballradius, boundary);
    
    %Setting those nodes to be at the lower plate (not below it):
    U(3*nodes,i+1) = uz_low;
    %Saying their velocity is now zero:
    v(3*nodes) = 0;
    
    %Updating force:
    f_last = f_current;
end


