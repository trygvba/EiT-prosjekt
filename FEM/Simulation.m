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

%[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);
[A M] = HomogenousMaterial(tetr(:,1:4),p,Cp, rhop);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V D] = eigs(A,M,20,'sm');
eigenvalues = diag(D);
%Pick the n'th eigenvalue we should analyse:
n = 8;
%omega = eigenvalues(n);
u = V(:,n);

alpha=0.5;


%Newmark 2beta method
dt=10^-9;
t_span=3*10^-8;
steps=100;     %Number of time increments.
sz=size(u,1);
              %Temporal step spacing.
U=zeros(sz,steps);
beta=0.25;
NumberOfPics = 100; %BEWARE OF NUMBER OF VTK-FILES.


K1=-M\A;
K2=(eye(sz) - beta*dt^2*K1)\eye(sz); 

%Initial data

v=zeros(sz,1);
U(:,1)=alpha*u;

output_folder = 'paraview/animation';
title = 'testing';

%Dirichlet boundary conditions

D=0.01; %amplitude of displacement
omega=2*pi/10; %frequency of displacement

dp =@(D,t,omega)  1.02- D*cos(omega*t);





n=1;
for i = 1:(steps-1)
    if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
        State_to_vtk(output_folder,title,n,sz,tetr(:,1:4),p,U(:,i));
        n = n+1;
    end
    U(:,i+1)=K2*(U(:,i) + dt*v + ((1- 2*beta)/2)*dt^2*K1*U(:,i));
    v=v + (dt/2)*K1*(U(:,i) +U(:,i+1));
    x=i/steps;
    waitbar(x)
end



    
