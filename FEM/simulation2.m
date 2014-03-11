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
%[A M] = HomogenousMaterial(tetr(:,1:4),p,Cp);



% TIME INTEGRATION, FORWARD EULER

Tspan=1;
T0=0.1;
szU=size(A,1);
steps=2;
dt=Tspan/steps;
D=0.5;
ballradius=max(p(:,3));
omega=2*pi;
u0=zeros(szU,1);
U=zeros(szU,steps);
U(:,1)=u0;
howlow=-1.04;
dp =@(D,t,omega)  1.20- D*cos(omega*t);
dU=zeros(szU,1);

for i=2:steps
   
 t=T0+i*dt;   

 
 [lowernodes,uzil] = lowerdirichletnodes( p,U(:,i-1), howlow, boundary );
 [ uppernodes,uzip ] = upperdirichletnodes( dp(D,t,omega), p, U(:,i-1), boundary );

 [fnew Anew Mnew unew dUnew] = IncorporateDirichletBoundary(A,M,U(:,i-1),dU,uppernodes,lowernodes,uzip,uzil);

 utemp= (eye(length(unew)) -(dt^2)/2*Mnew)\(fnew-Anew*unew) +dt*dUnew;
 U(:,i) = putDirichletBack(utemp, lowernodes, uppernodes, uzil, uzip);
 dU=(1/dt)*(U(:,i)-U(:,i-1));
 
end









% 
% 
% 
% 
% % vtk writing
% 
% output_folder = 'paraview/animation';
% title = 'testing';
% NumberOfPics = 100; %BEWARE OF NUMBER OF VTK-FILES.
% 
% n=1;
% for i = 1:(steps-1)
%     if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
%         State_to_vtk(output_folder,title,n,sz,tetr(:,1:4),p,U(:,i));
%         n = n+1;
%     end
%     U(:,i+1)=K2*(U(:,i) + dt*v + ((1- 2*beta)/2)*dt^2*K1*U(:,i));
%     v=v + (dt/2)*K1*(U(:,i) +U(:,i+1));
% end
% 
% %ODE45 versucht:
% % sz =size(u,1);
% % K1 = M\A;
% % F = @(t,u) [zeros(sz) K1; eye(sz) zeros(sz)]*u;
% % u0 =[zeros(sz,1); alpha*u];
% % tspan = [linspace(0,1,100)];
% % [TOUT YOUT] = ode45(F,tspan,u0);
% % U = YOUT(:,(sz+1):end)';
% 
% 
% 
% 
%     
