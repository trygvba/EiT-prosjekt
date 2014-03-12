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

Tspan=.001;
T0=0;
szU=size(A,1);
steps=20;
dt=Tspan/steps;
OLT=0.05;
impactzone=1;
ballradius=max(p(:,3));
omega=2*pi;

howlow=-1.04;
dp =@(OLT,t,omega)-OLT*sin(omega*t);
d2dp=@(OLT,t) OLT*omega^2*sin(omega*t);
u=zeros(3324,1);
maxdisplacement=ballradius-impactzone*OLT;
U=zeros(szU,steps);
dU=U(:,1);
[upperNodes, upperPlate ] = upperdirichletnodes( maxdisplacement, p, U(:,1), boundary );
[lowerNodes, lowerPlate] = lowerdirichletnodes( p,U(:,1), howlow, boundary );



upN=ones(length(upperNodes),1);

U(upperNodes,1) = upperPlate;
U(lowerNodes,1) = lowerPlate;


for i=2:steps


    t=T0+i*dt;   


    [fnew Anew Mnew] = IncorporateDirichletBoundaryV2(A,M,upperNodes,lowerNodes,upperPlate,lowerPlate,d2dp(OLT,t));

    upperPlate= upperPlate+(upN*(dp(OLT,t,omega)-ballradius));

    U(:,i)= U(:,i-1) +((dt^2)/2)*Mnew\(fnew-Anew*U(:,i-1)) +dt*dU;
    dU=(1/dt)*(U(:,i)-U(:,i-1));
    
    
end


% vtk writing

output_folder = 'paraview/animation';
title = 'testing';
NumberOfPics = 20; %BEWARE OF NUMBER OF VTK-FILES.

n=1;
for i = 1:(steps-1)
    if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
        State_to_vtk(output_folder,title,n,szU,tetr(:,1:4),p,U(:,i));
        n = n+1;
    end

end







