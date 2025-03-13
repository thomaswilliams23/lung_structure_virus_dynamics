


%Makes Figures 6(e) and (f) from the manuscript



%% setup

%set up path
addpath ../helper_functions
addpath plot_helpers


%data source
data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_lineages/gen_';


%loop specs
all_gens = [1,5];
rep_ind_to_plot = 1;

%set up parameters for plot
num_lineages = 4;
plotting_final_grid = 1;
generation_lengths = [100, 100, 100, 100, 100];

%set up cell IDs
target_id = 1;
eclipse_ids = 2:(2+num_lineages-1);
infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_lineages - 1);
infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_lineages-1);
infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_lineages - 1);
infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_lineages-1);
num_cell_types = infected_donor_dead_ids(end);


%save?
save_stuff = 0;



%% main loop


%loop over source positions (generation 1 or 5, i.e. left or right edge)
for gen_val = all_gens

    %load in data and plot
    load(strcat(data_source, num2str(gen_val), '/sim_data_rep_', num2str(rep_ind_to_plot), '.mat'))
    cell_grid = abs(sim_data_this_rep.final_grid_state);
    plot_grid; 
    
    %optionally save
    if save_stuff
        pause(1)
        exportgraphics(gcf, strcat('/branching_lineages_snapshot_gen_', num2str(gen_val), '.png'), 'Resolution', 300);
    end

end
