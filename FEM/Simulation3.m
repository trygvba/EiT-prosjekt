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
sz = size(A,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------PARAMETERS FOR TIME INTEGRATION:---------
%Radius of ball:
radius = max(p(:,3));
%Starting time:
T0 = 0;
%Time steps:
steps = 100;
%Temporal step size:
dt = 0.01;

%Upper Plate:
D = 0.1;
omega = 2*pi;
eps = 0.01;
plate = @(t) radius-eps-D*sin(omega*t);

%Lower Plate:
howlow = -(radius-eps);

%Initializings:
U = zeros(sz,steps);
v = sparse(sz,1);

[upnodes uz_up] = upperdirichletnodes(plate(T0),p,U(:,1),boundary);
[lownodes uz_low] = lowerdirichletnodes(p,U(:,1),howlow,boundary);

[ftemp Atemp Mtemp utemp vtemp] = IncorporateDirichletBoundary(A,M,U(:,1),v,upnodes,lownodes, uz_up, uz_low);
Utemp = putDirichletBack(utemp,lownodes,upnodes,uz_low,uz_up);