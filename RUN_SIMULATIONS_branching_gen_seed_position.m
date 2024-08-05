


%Conducts a sweep of simulations on the branching tree with infection 
%seeded by a single cell at a range of depths. Generates data of the type 
%used in Figures 2 and 3 of the manuscript.



%% setup

addpath helper_functions

%output folder
base_folder = 'SIMULATION_OUTPUT_branching_tree_seed_generation_sweep';
mkdir(base_folder)


%loop setup
num_reps = 10;
seed_depths = 1:25:500;

%parameter defaults
biology_parameters;
control_parameters;
prms = make_default_param_struct();


%specify parameters for the simulation
prms.vis_grids = 0;
prms.export_frames = 0;

prms.total_tissue_size_x = 5*100;
prms.total_tissue_size_y = 64;
prms.generation_lengths = [100, 100, 100, 100, 100];

prms.final_time = 300;
prms.stop_when_all_infected = 0;
prms.stop_when_all_dead = 1;



%specify data output times
vis_dt = 0.1;
out_times = 0:vis_dt:(prms.final_time-vis_dt);
num_out_times = length(out_times);



%% main loop
for seed_ind = 1:length(seed_depths)

    %specify x index of seed cell
    x_pos = seed_depths(seed_ind);

    %set up output data structures
    sim_data_this_seed_position{num_reps} = [];
    y_seed_positions = zeros(num_reps,1);


    %loop over iterations
    parfor rep = 1:num_reps
    
        %draw a random seed cell
        init_cell = [x_pos, randi(prms.total_tissue_size_y)];
    
        %run
        SIM_OUT = single_run_as_func(init_cell, out_times, prms);
    
        %save into arrays/structures
        y_seed_positions(rep) = init_cell(2);
        sim_data_this_seed_position{rep} = SIM_OUT;
    
    end

    %target folder setup
    target_folder = strcat(base_folder, '/x_seed_position_', num2str(x_pos));
    mkdir(target_folder)

    %save data
    save(strcat(target_folder, '/sim_data_this_seed_position'), 'sim_data_this_seed_position')
    save(strcat(target_folder, '/y_seed_positions'), 'y_seed_positions')

end


