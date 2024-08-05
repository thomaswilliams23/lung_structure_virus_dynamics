


% update extracellular virus



%% intialise

%check if needs ghost nodes based on boundary conditions
if ~use_toroidal_BCs
    new_virus=[v_ext(1,:,:); v_ext; v_ext(end,:,:)];
else 
    new_virus = v_ext;
end



%% implicit diffusion

%solve
for species = 1:num_lineages
    old_v_as_vec = reshape(new_virus(: , : ,species),[total_mesh_nodes,1]);
    next_v_as_vec = diff_matrix\old_v_as_vec;
    new_virus(:, :, species) = reshape(next_v_as_vec, [mesh_width, mesh_length]);
end


%trim ghost nodes
if ~use_toroidal_BCs
    new_virus=new_virus(2:(end-1),:,:);
end


%% explicit reaction


%decay
new_virus = new_virus - dt*c*v_ext;


%production
production_ids = [infected_donor_ids; infected_recip_ids];
is_producing = 0*v_ext;
for species = 1:num_lineages
    for cell_type = 1:size(production_ids, 1)
        is_producing(:,:,species) = is_producing(:,:,species) + (cell_grid==production_ids(cell_type,species));
    end
end
new_virus= new_virus + (dt/cell_area) * p * (1/numel(cell_grid)) * is_producing;



%update
v_ext=new_virus;


