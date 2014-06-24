function [ Ux Uy Uz ] = group_nodes( ux,uy,uz,tetr,phys_group)

    i=phys_group;
   
    index=unique(tetr(find(tetr(:,5)==i),(1:4)));
    Ux=zeros(length(index),length(ux(1,:)));
    Uy=zeros(length(index),length(ux(1,:)));
    Uz=zeros(length(index),length(ux(1,:)));
    Ux=ux(index,:);
    Uy=uy(index,:);
    Uz=uz(index,:);
  
    
end

