function [ e E E_analytic ] = energy_norm( tetr,ux,uy,uz,phys_groups,dt,steps,OLT,omega,p,wvel)

e=zeros(3,length(phys_groups));

for j=1:length(phys_groups)
    
    [Ux,Uy,Uz]=group_nodes(ux,uy,uz,tetr,phys_groups(j));
    
    for k=1:length(Ux(:,1))
        

        [Fx fx]=discretefourier(Ux,steps,dt,k);
    
    
        e(1,j)=e(1,j)+sum((Fx.^2).*(fx.^2));
    
        [Fy fy]=discretefourier(Uy,steps,dt,k);
    
        e(2,j)=e(2,j)+sum((Fy.^2).*(fy.^2));
    
        [Fz fz]=discretefourier(Uz,steps,dt,k);
    
        e(3,j)=e(3,j)+sum((Fz.^2).*(fz.^2));
        
        
    end
end

%total energy
N=sum(p(:,3)==max(p(:,3)));

E=sum(sum(e))/length(p(:,1));

E_analytic=OLT^2*pi*omega^2*wvel;


end
