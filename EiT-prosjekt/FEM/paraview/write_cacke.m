function write_cacke(output_file_name, title, node_num, element_num, element_order, element_node, cord, temp)

  output_file_name
  output_file = fopen(output_file_name, 'w');  

  fprintf ( output_file, '# vtk DataFile Version 2.0\n' );
  fprintf ( output_file, '%s\n', title );
  fprintf ( output_file, 'ASCII\n' );
  fprintf ( output_file, '\n' );

  % Write out the grid
  fprintf ( output_file, 'DATASET UNSTRUCTURED_GRID\n' );
  fprintf ( output_file, 'POINTS %d double\n', node_num );

  for node = 1 : node_num
    fprintf ( output_file, '  %f  %f  %f\n', cord(node,:) );
  end

  
% Write out the cells
  cell_size = element_num * ( element_order + 1 );

  fprintf ( output_file, '\n' );
  fprintf ( output_file, 'CELLS  %d  %d\n', element_num, cell_size );
  
  for element = 1 : element_num
    fprintf ( output_file, '  %d', element_order );
    for order = 1 : element_order
      fprintf ( output_file, '  %d', element_node(element,order) - 1 );
    end
    fprintf ( output_file, '\n' );
  end
%
%
  fprintf ( output_file, '\n' );
  fprintf ( output_file, 'CELL_TYPES %d\n', element_num );

  for element = 1 : element_num
    fprintf ( output_file, '10\n' );
  end
  
  % Write out datapoints 
  fprintf ( output_file, '\n' );
  fprintf ( output_file, 'POINT_DATA %d\n', node_num );
  fprintf ( output_file, 'SCALARS temp double\n' );
  fprintf ( output_file, 'LOOKUP_TABLE default\n' );
  
  for node = 1 : node_num
    fprintf ( output_file, '  %f\n', temp(node) );
  end 

  % Write cell types 
  fprintf ( output_file, '\n' );
  fprintf ( output_file, 'CELL_DATA %d\n', element_num );
  fprintf ( output_file, 'SCALARS material INT\n' );
  fprintf ( output_file, 'LOOKUP_TABLE default\n' );
  
  for node = 1 : element_num
    fprintf ( output_file, '  %d\n', element_node(node,5) );
  end 
  
  fclose(output_file); 
