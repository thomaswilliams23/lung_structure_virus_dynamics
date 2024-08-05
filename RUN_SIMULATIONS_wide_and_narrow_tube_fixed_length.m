


%Conducts simulations on tubes of fixed length 500 cells and circumference
%4 cells or 64 cells. Generates data of the type used in Figure 8 of the
%manuscript.




%% setup

addpath helper_functions

%output folder
target_folder = 'SWEEP_OUTPUT_wide_and_narrow_tube_fixed_length';
mkdir(target_folder)


%loop setup
num_reps = 100;
tube_circ_sweep = [4,64];


%parameter defaults
biology_parameters;
control_parameters;
prms = make_default_param_struct();


%specify parameters for the simulation
prms.vis_grid = 0;
prms.export_frames = 0;

prms.final_time = 300;
prms.stop_when_all_infected = 0;

prms.total_tissue_size_x = 500;
prms.generation_lengths = total_tissue_size_x;



%specify data output times
vis_dt = 0.1;
out_times = 0:vis_dt:(prms.final_time-vis_dt);
num_out_times = length(out_times);




%% main loop
for tube_circ = tube_circ_sweep


    %grid setup
    prms.total_tissue_size_y = tube_circ;
    prms.total_cells = tube_circ*500;


    %output structures
    all_sim_data_this_fldr = {};
    all_sim_data_this_fldr{num_reps} = [];


    %output folder
    output_folder = strcat(target_folder, '/tube_circ_',num2str(tube_circ));
    mkdir(output_folder);


    %parallel loop
    parfor rep = 1:num_reps

        %initial cells
        init_cells = [1, randi(prms.total_tissue_size_y)];

        %run simulation
        SIM_OUT = single_run_as_func(init_cells, out_times, prms);

        %save output
        all_sim_data_this_fldr{rep} = SIM_OUT;

    end


    %save output
    save(strcat(output_folder, '/all_sim_data_this_fldr.mat'), 'all_sim_data_this_fldr');

end
