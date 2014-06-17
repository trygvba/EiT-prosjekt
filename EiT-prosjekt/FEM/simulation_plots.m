clear all;
close all;
format long;
addpath(genpath('../Converters'));


%Declaration of parameters: 
%Polymer: 
X = 1; %Length scale.

Ep = 2*10^6;  
vp = 0.3; 

rhop = 1.02*10; 


%Silver: 
Es = 80*10^9; 
vs = 0.35; 
rhos = 10.5;


%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Ep,vp,vp);

[p tri tetr] = loadGeo('spherewshell_thick');

boundary = unique(tri);
[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhop);


%% Time Integration

%--------------------------------------------
%Parameters for time integration:

T0=0;                       %start time.
szU=size(A,1);              %dimension of our system.
szP=szU/3;
                
%Steps etc

dt=1/(1*10^3);
t_max=500*dt;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       300*dt;
steps=ceil(t_max/dt)




OLT=0.05;                   %Outer Layer Thickness.
impactzone=.05;              %Parameter to decide which nodes are in the Dirichlet boundary.
ballradius=max(p(:,3));     %Total radius of the ball with outher shell.
lower_ballradius=min(p(:,3));
howlow=-ballradius;
beta=1/4;
f=0;
epsilon=0.05;









topdisp=[];
    

max_disp=0.01;
    
maxf=0;
MAXF=(2/3)*(Ep*sqrt(ballradius)/(1-vp^2))*(max_disp^(3/2));





minf=0;
period=t_max/9;
loadrate=MAXF/period;


%Setting initial velocity:
v = zeros(szU,1);
K1=-M\A;
K2=(eye(szU) - beta*dt^2*K1)\eye(szU); 




v=zeros(szU,1);
U = zeros(szU,steps);
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
h = waitbar(0,'Time integrating...')
for i=1:(steps-1)    
%     if i==1||(floor(Number_of_pics*i/steps)>floor(Number_of_pics*(i-1)/steps))
%         State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,i));
%         n = n+1;
%     end
      t = T0+i*dt;
     
      p_new=pupdate(p,U(:,i));
      U(:,i+1) = K2*(U(:,i)+dt*v+0.25*dt^2*K1*U(:,i)+0.25*dt^2*(M\(f_vec2(p_new,tri,plateForceValidering2(f,loadrate,maxf,minf,t-dt,period),epsilon,FT) +f_vec2(p_new,tri,plateForceValidering2(f,loadrate,maxf,minf,t,period),epsilon,FT))));
      v = v+0.5*dt*(M\(f_vec2(p_new,tri,plateForceValidering2(f,loadrate,maxf,minf,t-dt,period),epsilon,FT) +f_vec2(p_new,tri,plateForceValidering2(f,loadrate,maxf,minf,t,period),epsilon,FT)))+0.5*dt*K1*(U(:,i+1)+U(:,i));
      
      pfplot(i)=plateForceValidering2(f,loadrate,maxf,minf,t-dt,period);

      [nodes,uzi] = lowerdirichletnodes(p,U(:,i+1), howlow, boundary );
      U(3*nodes,i+1)=uzi;
      v(3*nodes)=0;
      nodes;

     x=i/steps;
     waitbar(x)
end
toc


index_top = find(p(:,3)==ballradius);
index_bottom=find(p(:,3)==lower_ballradius)
figure
plot(X*U(3*index_top,:));
figure
plot(X*U(3*index_bottom,:));


% output_folder = 'paraview/animation';
% title = 'testing';

figure
plot(pfplot)
 
%  for n=1:steps
%      State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,n));
%  end
