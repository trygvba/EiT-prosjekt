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

alpha=5;


%Newmark 2beta method
steps=100;
[sz kuk]=size(u);
dt=0.001;
U=zeros(sz,steps);
beta=0.25;

K1=M/A;
K2=(eye(sz) - beta*dt^2*K1)\eye(sz); 

%Initial data

v=zeros(sz,1);
U(:,1)=alpha*u;

for i = 1:(steps-1)
    U(:,i+1)=K2*(U(:,i) + dt*v + ((1- 2*beta)/2)*dt^2*K1*U(:,i));
    v=v + (dt/2)*K1*(U(:,i) +U(:,i+1));
end






output_folder = 'paraview/animation';
title = 'testing';



%Converting to prefered data structure:
for j=1:steps
    for i = 1:3:sz
        uvec(ceil(i/3),:) = [U(i,j) U(i+1,j) U(i+2,j)];
        mag(ceil(i/3),j) = (U(i,j)^2+U(i+1,j)^2+U(i+2,j)^2)^0.5;
    end
    %Writing the shite to VTK:
    %write_to_vtk(sprintf('%s/%s%d.vtk',output_folder,title,j),title, node_num, element_num, element_order, element_node, cord, temp(:,i))
    writeVTK([output_folder '/' title '_' num2str(j)],tetr(:,1:4),p+uvec,mag(:,j));
end



    
