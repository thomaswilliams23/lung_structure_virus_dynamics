


%set up parameters for the PDE



%cell geometry
cell_width = 1;
cell_side_length = (1/sqrt(3)) * cell_width;
cell_area = (3*sqrt(3)/2) * cell_side_length^2;


%mesh elements
dx = 1;


%mesh size (total)
mesh_width = round((total_tissue_size_x)/dx)+2*(~use_toroidal_BCs);
mesh_length = round(total_tissue_size_y/dx);
total_mesh_nodes = mesh_width*mesh_length;


%mesh generation sizes
mesh_generation_lengths = generation_lengths;
if ~use_toroidal_BCs
    mesh_generation_lengths(1) = mesh_generation_lengths(1)+1;
    mesh_generation_lengths(end) = mesh_generation_lengths(end)+1;
end

mesh_generation_widths = generation_widths;
