


%Makes Figures 5(e), (f), (g), (h) and 6(c) and (d) from the manuscript



%% setup

%set up path
addpath ../helper_functions
addpath plot_helpers


%data sources
tube_data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_tube_lineages_full_data/tube_circ_';
toroid_data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_toroid_lineages_full_data/tube_circ_';


%loop specs
all_circs = [4, 8, 16, 32, 64, 128];
rep_ind_to_plot = 1;

%set up parameters for plot
num_lineages = 4;
plotting_final_grid = 1;

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

%set up output folder   
if save_stuff
    target_folder = 'lineages_snapshots';
    mkdir(target_folder);
end


%loop over circumferences
for circ_val = all_circs

    %load in tube data and plot
    load(strcat(tube_data_source, num2str(circ_val), '/sim_data_', num2str(rep_ind_to_plot), '.mat'))
    cell_grid = abs(sim_data_this_rep.final_grid_state);
    generation_lengths = size(cell_grid, 1);
    plot_grid; 
    
    %optionally save
    if save_stuff
        pause(1)
        exportgraphics(gcf, strcat(target_folder, '/tube_snapshot_circ_', num2str(circ_val), '.png'), 'Resolution', 300);
    end



    %load in toroid data and plot
    load(strcat(toroid_data_source, num2str(circ_val), '/sim_data_', num2str(rep_ind_to_plot), '.mat'))
    cell_grid = abs(sim_data_this_rep.final_grid_state);
    generation_lengths = size(cell_grid, 1);
    plot_grid; 
    
    %optionally save
    if save_stuff
        pause(1)
        exportgraphics(gcf, strcat(target_folder, '/toroid_snapshot_circ_', num2str(circ_val), '.png'), 'Resolution', 300);
    end
end
