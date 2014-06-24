function [ n_d ]=node_density(phys_group,tetr,p,rho)

i=phys_group;

indexT=tetr(find(tetr(:,5)==i),(1:4));
indexU=unique(tetr(find(tetr(:,5)==i),(1:4)));

n_d = zeros(length(indexU),1);
%length(n_d)

for j=length(n_d)
    
   [T y]=find(indexU(j)==indexT);
   len = length(T)
   for k=1:length(T)
       
       n_d(j)=n_d(j)+tetmass2(p,tetr,rho, T(k));
   end
   n_d(j)=(1/4)*n_d(j)
   
end

        

end



