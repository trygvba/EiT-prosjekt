clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters: 
%Polymer: 
X = 1;%15*10^(-6); %Length scale. ;
Ep = 10^6% 0.8*10^9*X;  
vp = 0.46; 
rhop = 950*X^3; 

%Silver: 
Es = 72.4*10^9*X; 
vs = 0.37; 
rhos = 10^4*X^3;


%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Ep,vp,vp);

[p tri tetr] = loadGeo('spherewshell');
boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhop);

%% Time Integration

%--------------------------------------------
%Parameters for time integration:

T0=0;                       %start time.
szU=size(A,1);              %dimension of our system.
szP=szU/3;
steps=400;                  %Number of time steps.
U = zeros(szU,steps);
dt=1/(steps);               %Temporal step size.
OLT=0.05;                   %Outer Layer Thickness.
impactzone=.05;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=2*pi;                 %Frequency of upperplate.
howlow=-ballradius*0.99;
beta=1/4;
f=1000000;

disp('kym spiser lorde-suppe')

%Setting initial velocity:
v = zeros(szU,1);

%------------------------------------------------

plateDisp = @(t) ballradius-OLT*sin(omega*t);
plateAcc = @(t) OLT*omega^2*sin(omega*t);


%plateForce=@(t)  1000*omega^2*sin(omega*t);


K1=-M\A;
K2=(eye(szU) - beta*dt^2*K1)\eye(szU); 




v=zeros(szU,1);

%------------------------------------------------
%   TIME INTEGRATION (NEWMARK)
%------------------------------------------------



tic
for i=2:steps
    
      t = T0+i*dt;
      [nodes,uzi] = lowerdirichletnodes(p,U(:,i-1), howlow, boundary );
      U(3*nodes,i)=uzi;
      p_new=pupdate(p,U(:,i-1),szU);
      U(:,i) = K2*(U(:,i-1)+dt*v+0.25*dt^2*K1*U(:,i-1)+0.25*dt^2*(M\(f_vec(p_new,tri,plateForce3(t-dt,omega,OLT,f),plateDisp(t-dt)) +f_vec(p_new,tri,plateForce3(t,omega,OLT,f),plateDisp(t)))));
      v = v+0.5*dt*(M\(f_vec(p_new,tri,plateForce3(t-dt,omega,OLT,f),plateDisp(t-dt)) +f_vec(p_new,tri,plateForce3(t,omega,OLT,f),plateDisp(t))))+0.5*dt*K1*(U(:,i-1)+U(:,i));
      
      
% %forward euler
% 
% U(:,i)=U(:,i-1) + (0.5*dt^2)*(M\(f_vec(p,tri,-plateForce(t),plateDisp(t)) - (A\U(:,i-1))))  + v*dt;
% v=(1/dt)*(U(:,i)-U(:,i-1));
%     
end
toc  


%------------POST PROCESSING:-----------------
output_folder = 'paraview/animation';
title = 'testing';

for n=1:steps
    State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,n));
end


