clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters: 
%Polymer: 
X = 15*10^(-6)/2; %Length scale.
Ep = 2*10^9*X;  
vp = 0.3; 

rhop = 1.02*10^3*X^3; 

%Silver: 
Es = 80*10^9*X; 
vs = 0.35; 
rhos = 10.5*10^3*X^3;


%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

[p tri tetr] = loadGeo('spherewshell');

boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V D] = eigs(A,M,20,'sm');
eigenvalues = sqrt(diag(D));
%Pick the n'th eigenvalue we should analyse:
n = 8;
%omega = eigenvalues(n);
u = V(:,n);
u = u/sqrt(u'*u);

alpha=0.1;
phase = @(t) alpha*cos(eigenvalues(n)*t);

%Newmark 2beta method
dt=10^-9;
t_span=3*10^-7;
steps=t_span/dt;     %Number of time increments.
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
Uana = U;
output_folder = 'paraview/animation';
title = 'testing';

err=zeros(1,steps-1);



n=1;
for i = 1:(steps-1)
    if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
        State_to_vtk(output_folder,title,n,sz,tetr(:,1:4),p,U(:,i));
        n = n+1;
    end
    U(:,i+1)=K2*(U(:,i) + dt*v + ((1- 2*beta)/2)*dt^2*K1*U(:,i));
    Uana(:,i+1) = Uana(:,1)*phase(i*dt);
    v=v + (dt/2)*K1*(U(:,i) +U(:,i+1));
    
    err(i)=max(abs(U(:,i+1)-Uana(:,i+1)));
    
    x=i/steps;
    waitbar(x)
end

plot(err)



