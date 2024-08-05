


%sets up the diffusion matrix as a sparse matrix



%% setup

%get i indices: one for each node in the system, seven times (self and six
%neighbours
run_up_nodes = 1:total_mesh_nodes;
i_inds = repmat(run_up_nodes, [1, 7]);

%values associated with each of these indices
diff_vals = [(1+4*virus_diff*dt/dx^2)*ones(1,total_mesh_nodes), ...
    -(2/3)*(virus_diff*dt/dx^2)*ones(1,total_mesh_nodes), ...
    -(2/3)*(virus_diff*dt/dx^2)*ones(1,total_mesh_nodes), ...
    -(2/3)*(virus_diff*dt/dx^2)*ones(1,total_mesh_nodes), ...
    -(2/3)*(virus_diff*dt/dx^2)*ones(1,total_mesh_nodes), ...
    -(2/3)*(virus_diff*dt/dx^2)*ones(1,total_mesh_nodes), ...
    -(2/3)*(virus_diff*dt/dx^2)*ones(1,total_mesh_nodes)];


%% find j_inds (i.e. neighbours of each cell in linear indices)

%initialise
bl_neighbour = zeros(1,total_mesh_nodes);
tl_neighbour = zeros(1,total_mesh_nodes);
b_neighbour = zeros(1,total_mesh_nodes);
t_neighbour = zeros(1,total_mesh_nodes);
br_neighbour = zeros(1,total_mesh_nodes);
tr_neighbour = zeros(1,total_mesh_nodes);

grid_dims = [mesh_width, mesh_length];



%% read neighbours in

%loop over nodes
for node_i = 1:mesh_width
    for node_j = 1:mesh_length

        %get linear index
        cell_ind = sub2ind(grid_dims, node_i, node_j);

        %get indices of neighbours
        [neighbour_i, neighbour_j] = get_neighbour_indices(...
            node_i, node_j, grid_dims,  mesh_generation_lengths, mesh_generation_widths);

        %get linear indices of neighbours
        bl_neighbour(cell_ind) = sub2ind(grid_dims, neighbour_i(1), neighbour_j(1));
        tl_neighbour(cell_ind) = sub2ind(grid_dims, neighbour_i(2), neighbour_j(2));
        b_neighbour(cell_ind) = sub2ind(grid_dims, neighbour_i(3), neighbour_j(3));
        t_neighbour(cell_ind) = sub2ind(grid_dims, neighbour_i(4), neighbour_j(4));
        br_neighbour(cell_ind) = sub2ind(grid_dims, neighbour_i(5), neighbour_j(5));
        tr_neighbour(cell_ind) = sub2ind(grid_dims, neighbour_i(6), neighbour_j(6));

    end
end


j_inds = [run_up_nodes, bl_neighbour, tl_neighbour, b_neighbour, t_neighbour, br_neighbour, tr_neighbour];



%% specify j_inds for left and right edge nodes (enforces toroidal geometry)
    
%left edge
node_i = 1;
for node_j = 1:mesh_length

    %get linear index
    cell_ind = sub2ind(grid_dims, node_i, node_j);

    %add in the relevant j indices
    j_inds(1*total_mesh_nodes + cell_ind) = sub2ind(grid_dims, mesh_length, mod(node_j-2,mesh_length)+1); %BL
    j_inds(2*total_mesh_nodes + cell_ind) = sub2ind(grid_dims, mesh_length, node_j); %TL
end


%right edge
node_i = mesh_width;
for node_j = 1:mesh_length

    %get linear index
    cell_ind = sub2ind(grid_dims, node_i, node_j);

    %add in the relevant j indices
    j_inds(5*total_mesh_nodes + cell_ind) = sub2ind(grid_dims, 1, node_j); %BR
    j_inds(6*total_mesh_nodes + cell_ind) = sub2ind(grid_dims, 1, mod(node_j, mesh_length)+1); %TR

end


%% no-flux case (override toroidal setup)
if ~use_toroidal_BCs

    %left edge (ghosts)
    node_i = 1;
    for node_j = 1:mesh_length
    
        %get linear index
        cell_ind = sub2ind(grid_dims, node_i, node_j);
    
        %fix coeffs:
        diff_vals(1*total_mesh_nodes + cell_ind) = 0; %BL
        diff_vals(2*total_mesh_nodes + cell_ind) = 0; %TL
    
        diff_vals(5*total_mesh_nodes + cell_ind) = 2*diff_vals(5*total_mesh_nodes + cell_ind);%BR
        diff_vals(6*total_mesh_nodes + cell_ind) = 2*diff_vals(6*total_mesh_nodes + cell_ind);%TR
    end
    
    
    %right edge (ghosts)
    node_i = mesh_width;
    for node_j = 1:mesh_length
    
        %get linear index
        cell_ind = sub2ind(grid_dims, node_i, node_j);
    
        %fix coeffs:
        diff_vals(5*total_mesh_nodes + cell_ind) = 0; %BR
        diff_vals(6*total_mesh_nodes + cell_ind) = 0; %TR
    
        diff_vals(1*total_mesh_nodes + cell_ind) = 2*diff_vals(1*total_mesh_nodes + cell_ind);%BL
        diff_vals(2*total_mesh_nodes + cell_ind) = 2*diff_vals(2*total_mesh_nodes + cell_ind);%TL
    
    end
end




%% set up diffusion matrix
sparse_size = length(diff_vals);
diff_matrix = sparse(i_inds, j_inds, diff_vals, total_mesh_nodes, total_mesh_nodes, sparse_size);
