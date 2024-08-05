


%Conducts a sweep of simulations on tubes of varying circumferences and
%with varying viral diffusion coefficients. Generates data of the type used
%in Figure 1 of the manuscript.



%% setup

addpath helper_functions

%output folder
target_folder = 'SWEEP_OUTPUT_diffusion_and_aspect_ratio_sweep';
mkdir(target_folder)


%loop setup
num_reps = 10;
tube_circ_sweep = 4; %[4,8,16,32,64];
diff_sweep = [0, 1]; %[0, 1, 10, 100, 1000, Inf];


%fixed number of cells
total_cells = 64^2;



%parameter defaults
biology_parameters;
control_parameters;
prms = make_default_param_struct();


%specify parameters for the simulation
prms.vis_grid = 0;
prms.export_frames = 0;

prms.final_time = 6000; %3000;
prms.stop_when_all_infected = 1;
prms.stop_if_dieout = 1;



%specify data output times
vis_dt = 0.1;
out_times = 0:vis_dt:(prms.final_time-vis_dt);
num_out_times = length(out_times);




%% main loop
for tube_circ = tube_circ_sweep

    %% loop diffusion
    for diff_val = diff_sweep


        %diffusion setup
        prms.virus_diff = diff_val;

        %grid setup
        prms.total_tissue_size_y = tube_circ;
        prms.total_tissue_size_x = round(total_cells/tube_circ);
        prms.generation_lengths = prms.total_tissue_size_x;
        prms.total_cells = total_cells;


        %output structures
        all_sim_data_this_fldr = {};
        all_sim_data_this_fldr{num_reps} = [];


        %output folder
        output_folder = strcat(target_folder, '/tube_circ_',num2str(tube_circ),...
            '_diff_', num2str(diff_val));
        mkdir(output_folder);


        %parallel loop
        parfor rep = 1:num_reps

            %initial cells
            num_init_infected = 4;
            init_cells = ones(num_init_infected, 2);
            init_cells(:,2) = randperm(prms.total_tissue_size_y,num_init_infected);

            %run simulation
            SIM_OUT = single_run_as_func(init_cells, out_times, prms);

            %save output
            all_sim_data_this_fldr{rep} = SIM_OUT;

        end


        %save output
        save(strcat(output_folder, '/all_sim_data_this_fldr.mat'), 'all_sim_data_this_fldr');

    end

end
