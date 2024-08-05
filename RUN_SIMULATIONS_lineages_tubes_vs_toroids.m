


%Conducts a sweep of simulations on tubes of varying circumferences with 
%fixed length, and on toroids with an approxmately equal number of cells 
%and equal aspect ratio. Generates data of the type used in Figure 4 of 
%the manuscript.



%% setup

addpath helper_functions

%output folders
tube_target_folder = 'SIMULATION_OUTPUT_tube_circ_sweep';
mkdir(tube_target_folder)
toroid_target_folder = 'SIMULATION_OUTPUT_toroid_circ_sweep';
mkdir(toroid_target_folder)


%loop setup
num_reps = 100;
tube_circ_sweep = [4,8,16,32,64,128];


%parameter defaults
biology_parameters;
control_parameters;
tube_prms = make_default_param_struct();


%specify parameters for the simulation
tube_prms.vis_grid = 0;
tube_prms.export_frames = 0;
tube_prms.final_time = 200;
tube_prms.num_lineages = 4;
tube_prms.total_tissue_size_x = 500;
tube_prms.generation_lengths = total_tissue_size_x;


%specify data output times
vis_dt = 0.1;
out_times = 0:vis_dt:(tube_prms.final_time-vis_dt);
num_out_times = length(out_times);




%% loop tube circumferences
for tube_circ = tube_circ_sweep


    %% first, the tube case

    %grid setup
    tube_prms.total_tissue_size_y = tube_circ;
    tube_prms.total_cells = tube_prms.total_tissue_size_x * tube_prms.total_tissue_size_y;

    %output structures
    all_sim_data = {};
    all_sim_data{num_reps} = [];


    %output folder
    output_folder = strcat(tube_target_folder, '/tube_circ_',num2str(tube_circ));
    mkdir(output_folder);


    %parallel loop
    for rep = 1:num_reps

        %initial cells
        num_init_infected = 4;
        init_cells = ones(num_init_infected, 2);
        init_cells(:,2) = randperm(tube_prms.total_tissue_size_y,num_init_infected);

        %run simulation
        SIM_OUT = single_run_as_func(init_cells, out_times, tube_prms);

        %save output
        all_sim_data{rep} = SIM_OUT;

    end


    %save output one sim at a time (to avoid having huge data structures)
    for rep = 1:num_reps
        sim_data_this_rep = all_sim_data{rep};
        save(strcat(output_folder,'/sim_data_',num2str(rep),'.mat'), 'sim_data_this_rep');
    end




    %% then, the toroid case

    %make a new parameter struct
    toroid_prms = tube_prms;

    %fix the grid setup
    toroid_prms.total_tissue_size_x = round(srqt(tube_prms.total_cells));
    toroid_prms.total_tissue_size_y = round(sqrt(tube_prms.total_cells));
    toroid_prms.total_cells = toroid_prms.total_tissue_size_x*toroid_prms.total_tissue_size_y;

    %impose toroidal boundaries
    toroid_prms.use_toroidal_BCs = 1;


    %output structures
    all_sim_data = {};
    all_sim_data{num_reps} = [];


    %output folder
    output_folder = strcat(toroid_target_folder, '/toroid_circ_',num2str(tube_circ));
    mkdir(output_folder);


    %parallel loop
    for rep = 1:num_reps

        %initial cells
        num_init_infected = 4;
        init_cells = ones(num_init_infected, 2);
        init_cells(:,1) = randperm(tube_prms.total_tissue_size_x,num_init_infected);
        init_cells(:,2) = randperm(tube_prms.total_tissue_size_y,num_init_infected);

        %run simulation
        SIM_OUT = single_run_as_func(init_cells, out_times, toroid_prms);

        %save output
        all_sim_data{rep} = SIM_OUT;

    end


    %save output one sim at a time (to avoid having huge data structures)
    for rep = 1:num_reps
        sim_data_this_rep = all_sim_data{rep};
        save(strcat(output_folder,'/sim_data_',num2str(rep),'.mat'), 'sim_data_this_rep');
    end

end
