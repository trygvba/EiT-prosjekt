function [ n_d ]=node_density(phys_group,tetr,p,rho)

i=phys_group;

indexT=tetr(find(tetr(:,5)==i),(1:4));
indexU=unique(indexT);
%plot(indexU,'*')
%size(indexU)

n_d = zeros(length(indexU),1);
%length(n_d)

for j=1:length(n_d)
    
   [T y]=find(indexT==indexU(j));
   %len = size([T y])
   for k=1:length(T)
       
       n_d(j)=n_d(j)+tetmass2(p,indexT,rho, T(k))/4;
   end
   
end

        

end



    