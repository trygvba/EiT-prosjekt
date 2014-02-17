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
Es = 72.4*10^9
vs = 0.37;
rhos = 10*10^3;


%-----------ASSEMBLY:------------------------
[Cp Cs] = StressMatrices(Ep,Es,vp,vs);

[p tri tetr] = loadGeo('spherewshell');

[A M] = MassAndStiffnessMatrix3D(tetr,p,Cp,Cs,rhop,rhos);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V D] = eigs(A,M,20,'sm');
eigenvalues = diag(D);
%Pick the n'th eigenvalue we should analyse:
n = 8;
omega = eigenvalues(n);
u = V(:,n);

alpha=1;

%Converting to prefered data structure:
for i = 1:3:length(u)
	uvec(ceil(i/3),:) = [u(i) u(i+1) u(i+2)];
    mag(ceil(i/3),:) = (u(i)^2+u(i+1)^2+u(i+2)^2)^0.5;
end

writeVTK('test',tetr(:,1:4),p,mag);
%writeVTK(['../PostProc/eigenmode_' num2str(n)],tetr,p,mag);
