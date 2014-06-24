%clear all;
%close all;
format long;
addpath(genpath('../Converters'));

%Declaration of parameters:
%Polymer:
X = 1;%15*10^(-6); %Length scale.
Phys_groups = [70,71,72];
E = 10^(0)*[1,1,1];%10^9*X;
v = [0,0,0];
rho = [4,.5,4];%X^(3)*950;

Meshname = 'QMembrane'; % I mesh-detail-rekkefølge: HQMemsphere,SMemSphere Membranewsphere

harmonicMode = 10; % Harmonic mode (2n+1)
ispulse = 1;

%Parameters for Paraview printing:
ExtraNameNote = 'omega_';

modestring = sprintf('%d',harmonicMode);
vtktitle = [Meshname '_' ExtraNameNote modestring(1:(min(3,length(modestring)))) ];
if ispulse
    vtktitle = [vtktitle 'p' num2str(ispulse)];
end
vtktitle

output_folder = 'paraview/animation/Membrane';
steps=25;    % Number of time steps.
NumberOfPics = steps;
granularity = 0.1;
OLT = 0.05;



%Parameters for time integration:
T0=0;                       %start time.
%szU=size(A,1);              %dimension of our system.
tic
[p, tri, tetr, szU, K1, K2, Amod, Mmod_inv, vel,Fmat_up, Fmat_low, F_acc, lowerNodes, uz_low,upperNodes, uz_up, omega, dt, x_plates, y_plates] = Assembly_Membrane(Meshname,Phys_groups, E,v,rho, harmonicMode, granularity);
toc
plateDisp = @(t) -OLT*((sin(omega*t))).*(ispulse*t<(2*pi/omega));
plateAcc = @(t) -omega^2*plateDisp(t);

wvel = omega/(harmonicMode*pi/4);
%AnalyticSolution = @(x,t) plateDisp(t+(x-1)/wvel) - plateDisp(t-(x+3)/wvel).*(wvel*t>=(x+3)) + plateDisp(t+(x-5)/wvel).*(wvel*t>=(5-x)) - plateDisp(t-(x+7)/wvel).*(wvel*t>=(7+x));

MarkerNode = find(sum(p.^2,2)==min(sum(p.^2,2)));

TimeIntegrator;
% F_last = -Fmat_up*uz_up-Fmat_low*uz_low;
% Utemp_last = Amod\F_last;
% u0 = putDirichletBackV2(Utemp_last,lowerNodes,upperNodes,uz_low,uz_up,x_plates,y_plates);
% %U = full(u0);
% U = zeros(szU,1);
% disp('Starting time integration')
% %------------------------------------------------
% %       TIME INTEGRATION (NEWMARK)
% %------------------------------------------------
% Uz=zeros(szU/3,steps);
% n=1;
% %h = waitbar(0,'Pictures taken');
% tic
% for i=1:steps
%     
% %     if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
% %         State_to_vtk(output_folder,vtktitle,n,szU,tetr(:,1:4),p,U);
% %         n=n+1;
% %     end
%     t =T0+(i-1)*dt;
%     F_cur = -Fmat_up*(uz_up+plateDisp(t))-Fmat_low*uz_low-plateAcc(t)*F_acc;
%     Utemp_cur = K2*(Utemp_last+dt*vel-0.25*dt^2*K1*Utemp_last+0.25*dt^2*(Mmod_inv*(F_cur+F_last)));
%     utemp = putDirichletBackV2(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t),x_plates,y_plates);
%     U = utemp;
%     Uz(:,i) = U(3:3:end);
%     vel = vel+ 0.5*dt*(Mmod_inv*(F_cur+F_last))-0.5*dt*K1*(Utemp_cur+Utemp_last);
%     Utemp_last = Utemp_cur;
%     F_last = F_cur;
%     %waitbar(i/steps); 
% end
% toc
%close(h);
disp('Done');