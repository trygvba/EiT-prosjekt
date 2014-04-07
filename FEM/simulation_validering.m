clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters: 
%Polymer: 
X = 15*10^(-6); %Length scale.
Ep = 3*10^9*X;  
vp = 0.48; 

rhop = 950*X^3; 

%Silver: 
Es = 83*10^9*X; 
vs = 0.37; 
rhos = 10490*X^3;


%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

[p tri tetr] = loadGeo('spherewshell');

boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);


%% Time Integration

%--------------------------------------------
%Parameters for time integration:

T0=0;                       %start time.
szU=size(A,1);              %dimension of our system.
szP=szU/3;
                
%Steps etc
dt=1/(1*10^9);
t_max=8*10^-7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       300*dt;
steps=ceil(t_max/dt)


U = zeros(szU,steps);
OLT=0.05;                   %Outer Layer Thickness.
impactzone=.05;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
omega=10^12*pi;                 %Frequency of upperplate.
howlow=-ballradius;
beta=1/4;
f=-(10^-1);%/X;
epsilon=0.01;
disp('kym spiser lorde-suppe')


maxf=10^-1/X;


minf=0;
period=t_max/5;
loadrate=-maxf/period;


%Setting initial velocity:
v = zeros(szU,1);
K1=-M\A;
K2=(eye(szU) - beta*dt^2*K1)\eye(szU); 




v=zeros(szU,1);

%------------------------------------------------
%   TIME INTEGRATION (NEWMARK)
%------------------------------------------------

FT=getNeumannBoundary(tri,p,ballradius)
ft=p(FT(1,:),:);

% output_folder = 'paraview/animation';
% title = 'testing';
% Number_of_pics = 200;
% n=1;
% 
pfplot=zeros(1,steps-1);
tic
for i=1:(steps-1)    
%     if i==1||(floor(Number_of_pics*i/steps)>floor(Number_of_pics*(i-1)/steps))
%         State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,i));
%         n = n+1;
%     end
      t = T0+i*dt;
     
      p_new=pupdate(p,U(:,i));
      U(:,i+1) = K2*(U(:,i)+dt*v+0.25*dt^2*K1*U(:,i)+0.25*dt^2*(M\(f_vec2(p_new,tri,plateForceValidering1(f,loadrate,maxf,minf,t-dt,period),epsilon,FT) +f_vec2(p_new,tri,plateForceValidering(f,loadrate,maxf,minf,t,period),epsilon,FT))));
      v = v+0.5*dt*(M\(f_vec2(p_new,tri,plateForceValidering1(f,loadrate,maxf,minf,t-dt,period),epsilon,FT) +f_vec2(p_new,tri,plateForceValidering1(f,loadrate,maxf,minf,t,period),epsilon,FT)))+0.5*dt*K1*(U(:,i+1)+U(:,i));
      
      pfplot(i)=plateForceValidering1(f,loadrate,maxf,minf,t-dt,period)*X;

      [nodes,uzi] = lowerdirichletnodes(p,U(:,i+1), howlow, boundary );
      U(3*nodes,i+1)=uzi;
      v(3*nodes)=0;
      nodes;

     x=i/steps;
     waitbar(x)
end
toc  
index_top = find(p(:,3)==ballradius);
figure
plot(X*U(3*index_top,:));
output_folder = 'paraview/animation';
title = 'testing';

figure
plot(pfplot)
% 
% for n=1:steps
%     State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,n));
% end

