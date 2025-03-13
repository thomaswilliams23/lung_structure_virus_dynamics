


%Makes Figure 2(c) from the manuscript (up to stochastic noise)



%% run a single iteration of the model

curr_dir =  cd('../');
addpath helper_functions/


%load default params
prms = make_default_param_struct;

%%% CHANGE PARAMETERS HERE

%export frames
prms.vis_grid = 0;
prms.export_frames = 1;

%number of lineages
prms.num_lineages = 1;

%timing
out_times = [0, 20, 40, 60];
prms.final_time = out_times(end);

%dimensions - single tube
prms.total_tissue_size_x = 128;
prms.total_tissue_size_y = 32;
prms.all_cells_total = prms.total_tissue_size_x * prms.total_tissue_size_y;
prms.generation_lengths = [128];
assert(sum(prms.generation_lengths) == prms.total_tissue_size_x);

%diffusion
prms.virus_diff = 100;
%%%


%initial cells
num_init_infected = 4;
init_cells = ones(num_init_infected, 2);
init_cells(:,2) = randperm(prms.total_tissue_size_y,num_init_infected);


%run simulation
[sim_out] = single_run_as_func(init_cells, out_times, prms);





%% plot frames from the model output
cd(curr_dir)
addpath ../helper_functions/

num_lineages = prms.num_lineages;
generation_lengths = prms.generation_lengths;

for frame_ind = 1:length(out_times)
    cell_grid = abs(squeeze(sim_out.grid_frames(:,:,frame_ind)));
    plot_grid;
    pause(1)
    exportgraphics(gcf, strcat('Fig1c_snapshot_', num2str(frame_ind), '.png'), "Resolution",300)
    close
end




