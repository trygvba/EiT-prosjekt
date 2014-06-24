function [ Ux Uz Uy  ] = untitled3( ux,uy,uz,tetr,phys_group)

    i=phys_group;
   
    index=unique(tetr(find(tetr(:,5)==i),(1:4)));
    Ux=zeros(length(index(:,1)),length(ux(1,:)));
    Uy=zeros(length(index(:,1)),length(ux(1,:)));
    Uz=zeros(length(index(:,1)),length(ux(1,:)));
    Ux(1,:)=ux(index,:);
    Uy(2,:)=uy(index,:);
    Uz(3,:)=uz(index,:);
  
    
end

