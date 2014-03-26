function State_to_vtk(output_folder,title,n,size,tetr,p,U)
    for i = 1:3:size
        uvec(ceil(i/3),:) = [U(i) U(i+1) U(i+2)];
        mag(ceil(i/3)) = (U(i)^2+U(i+1)^2+U(i+2)^2)^0.5;
    end
    %Writing the shite to VTK:
    %write_to_vtk(sprintf('%s/%s%d.vtk',output_folder,title,j),title, node_num, element_num, element_order, element_node, cord, temp(:,i))
    writeVTK([output_folder '/' title '_' num2str(n)],tetr,p+uvec,mag);
end