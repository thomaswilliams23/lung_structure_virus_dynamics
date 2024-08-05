function default_prm = make_default_param_struct()

    %uses params in the parameter files to populate a param struct
    
    biology_parameters;
    control_parameters;
    mesh_parameters;
    
    
    %read in
    default_prm.total_tissue_size_x = total_tissue_size_x;
    default_prm.total_tissue_size_y = total_tissue_size_y;
    default_prm.all_cells_total = all_cells_total;
    
    default_prm.generation_lengths = generation_lengths;
    
    default_prm.final_time = final_time;
    default_prm.stop_when_all_infected = stop_when_all_infected;
    default_prm.stop_if_dieout = stop_if_dieout;
    
    default_prm.num_lineages = num_lineages;

    default_prm.use_toroidal_BCs = use_toroidal_BCs;
    
    default_prm.vis_grid = vis_grid;
    default_prm.export_frames = export_frames;
    default_prm.frame_interval = frame_interval;


    default_prm.latent_stages = latent_stages;
    default_prm.virus_diff = virus_diff;
    default_prm.alpha_param = alpha_param; 
    default_prm.beta_param = beta_param;
    default_prm.gamma_param = gamma_param;
    default_prm.delta_param = delta_param;  
    default_prm.p = p;
    default_prm.c = c;
    
end