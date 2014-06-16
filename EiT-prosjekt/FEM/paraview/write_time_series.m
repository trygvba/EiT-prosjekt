function write_time_series(output_folder, title, node_num, element_num, element_order, element_node, cord, temp, timesteps)
filecount = 1;
for i=1:timesteps:length(temp)
    write_to_vtk(sprintf('%s/%s%d.vtk',output_folder,title,filecount),title, node_num, element_num, element_order, element_node, cord, temp(:,i))
    filecount = filecount +1;
end
