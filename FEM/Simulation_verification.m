%clear all;
%close all;
format long;
addpath(genpath('../Converters'));

%Declaration of parameters:
%Polymer:
X = 1;%15*10^(-6); %Length scale.
Ep = 10^(5);%10^9*X;
vp = 0.49;
rhop = 1;%X^(3)*950;
vtktitle = 'PeriodMode5';
harmonicMode = 5; % Harmonic mode (2n+1)
steps=1000;        % Number of time steps.
OLT = 0.05;
MarkerNode = 562;

%Parameters for time integration:
T0=0;                       %start time.
%szU=size(A,1);              %dimension of our system.
[p, tri, tetr, szU, K1, K2, Amod, Mmod_inv, v,Fmat_up, Fmat_low, F_acc, lowerNodes, uz_low,upperNodes, uz_up, omega, dt, x_plates, y_plates] = Assembly_verification('QCube',Ep,vp, rhop, harmonicMode);
plateDisp = @(t) -OLT*sin(omega*t).*(t<(2*pi/omega));
plateAcc = @(t) OLT*omega^2*sin(omega*t);
wvel = omega/(0.96*harmonicMode*pi/4);
%AnalyticSolution = @(x,t) plateDisp(t+(x-1)/wvel) - plateDisp(t-(x+3)/wvel).*(wvel*t>=(x+3)) + plateDisp(t+(x-5)/wvel).*(wvel*t>=(5-x)) - plateDisp(t-(x+7)/wvel).*(wvel*t>=(7+x));

%))
%Parameters for Paraview printing:
output_folder = 'paraview/animation/QCube';
NumberOfPics = 10;

F_last = -Fmat_up*uz_up-Fmat_low*uz_low;
Utemp_last = Amod\F_last;
u0 = putDirichletBackV2(Utemp_last,lowerNodes,upperNodes,uz_low,uz_up,x_plates,y_plates);
%U = full(u0);
U = zeros(szU,1);
disp('Starting time integration')
%------------------------------------------------
%       TIME INTEGRATION (NEWMARK)
%------------------------------------------------
Uz=zeros(szU/3,steps);
AUz = Uz;
DUz = Uz;
n=1;
h = waitbar(0,'Pictures taken');
tic
for i=1:steps
    
    if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
        State_to_vtk(output_folder,vtktitle,n,szU,tetr(:,1:4),p,U);
        waitbar(n/NumberOfPics); 
        n = n+1;
    end
    t =T0+(i-1)*dt;
    F_cur = -Fmat_up*(uz_up+plateDisp(t))-Fmat_low*uz_low-plateAcc(t)*F_acc;
    Utemp_cur = K2*(Utemp_last+dt*v-0.25*dt^2*K1*Utemp_last+0.25*dt^2*(Mmod_inv*(F_cur+F_last)));
    utemp = putDirichletBackV2(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t),x_plates,y_plates);
    U = utemp;
    Uz(:,i) = U(3:3:end);
    AUz(:,i) = AnalyticSolver(p(:,3),(i-1)*dt, omega, wvel);
    v = v+ 0.5*dt*(Mmod_inv*(F_cur+F_last))-0.5*dt*K1*(Utemp_cur+Utemp_last);
    Utemp_last = Utemp_cur;
    F_last = F_cur;
end
toc
close(h);
hold on
plot(1:steps, Uz(MarkerNode,:))
plot(1:steps, AUz(MarkerNode,:),'r')
legend('Simulation', 'Analytic');
title('Displacement of node 562');
xlabel('Timesteps');
ylabel('Displacement relative original position.');
disp('Done');