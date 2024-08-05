


%Runs simulations on the branching tree seeded on either the left (open) 
%or right (branched) edge. Generates data of the type used in Figures 5-8 
%of the manuscript.



%% setup

addpath helper_functions



%loop setup
num_reps = 100;

%parameter defaults
biology_parameters;
control_parameters;
prms = make_default_param_struct();


%specify parameters for the simulation
prms.vis_grids = 0;

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



%%%%% THESE OPTIONS DEPEND ON APPLICATION %%%%%

%this version for lineage simulations (Figure 5-6)
prms.num_lineages = 4;
prms.export_frames = 1;
base_folder = 'SIMULATION_OUTPUT_branching_tree_lineages_left_and_right_edge'; 

%this version for applying immune action (Figure 7-8)
% prms.num_lineages = 1;
% prms.export_frames = 0;
% base_folder = 'SIMULATION_OUTPUT_branching_tree_immune_left_and_right_edge';

mkdir(base_folder)

%%%%%                                    %%%%%



%% main loop
for start_on_left_edge = [1,0]

    %specify x index of seed cell
    if start_on_left_edge
        x_pos = 1;
    else
        x_pos = prms.total_tissue_size_x;
    end

    %set up output data structures
    all_sim_data{num_reps} = [];


    %loop over iterations
    parfor rep = 1:num_reps
    
        %initial cells
        num_init_infected = prms.num_lineages;

        if num_init_infected==1

            %draw one random seed cell
            init_cells = [x_pos, randi(prms.total_tissue_size_y)];
        else

            %four initial cells (taking care they are in the same branch)
            init_cells = x_pos*ones(num_init_infected, 2);

            num_final_tubes = 2^(length(prms.generation_lengths)-1);
            width_of_final_tube = prms.total_tissue_size_y/num_final_tubes;
            start_of_rand_tube = (randi(num_final_tubes)-1)*width_of_final_tube;

            for cell_ind = 1:prms.num_lineages
                init_cells(cell_ind,2) = start_of_rand_tube + cell_ind;
            end
        end

    
        %run
        SIM_OUT = single_run_as_func(init_cells, out_times, prms);
    
        %save into arrays/structures
        all_sim_data{rep} = SIM_OUT;
    
    end

    
    %target folder setup
    if start_on_left_edge
        target_folder = strcat(base_folder, '/gen_1', num2str(x_pos));
    else
        target_folder = strcat(base_folder, '/gen_5', num2str(x_pos));
    end
    mkdir(target_folder)


    %save output one sim at a time (to avoid having huge data structures)
    for rep = 1:num_reps
        sim_data_this_rep = all_sim_data{rep};
        save(strcat(target_folder,'/sim_data_rep_',num2str(rep),'.mat'), 'sim_data_this_rep');
    end

end


