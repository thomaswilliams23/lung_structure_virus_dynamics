


%Makes Figure 7(c) and (d) from the manuscript



%% setup

%first set text intepreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath ../helper_functions
addpath plot_helpers


%data sources
branching_data = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_lineages/gen_1';
tube_data = '../SIM_DATA_from_ms/SWEEP_OUTPUT_tube_lineages_full_data/tube_circ_64';


%parameters and specs
num_reps = 100;
num_lineages = 4;
generation_lengths = [100, 100, 100, 100, 100];
generation_circs = [64, 32, 16, 8, 4];
grid_size_x = sum(generation_lengths);
grid_size_y = generation_circs(1);


%specify band width
band_width = 10;

%save?
save_stuff = 1;



%% loop

%compute number of bands and initialise data structures
num_bands = grid_size_x - 2 - band_width + 1;
num_spec_in_band_branching = zeros(num_reps,num_bands);
num_spec_in_band_tube = zeros(num_reps,num_bands);
kappa_in_band_branching = zeros(num_reps,num_bands);
kappa_in_band_tube = zeros(num_reps,num_bands);


%loop over reps
for rep = 1:num_reps

    
    %load in the grid state
    load(strcat(branching_data, '/sim_data_rep_', num2str(rep), '.mat'));
    final_grid_branching = sim_data_this_rep.final_grid_state;

    load(strcat(tube_data, '/sim_data_', num2str(rep), '.mat'));
    final_grid_tube = sim_data_this_rep.final_grid_state;

    
    
    %annotate grid lineages
    get_lineage_num = @(id) get_infected_lineage_incl_dead(id, num_lineages);
    final_grid_lineages_branching = arrayfun(get_lineage_num, abs(final_grid_branching));
    final_grid_lineages_tube = arrayfun(get_lineage_num, abs(final_grid_tube));
    

    %loop over bands
    parfor band_ind = 1:num_bands
        

        %count number of species in band
        num_spec_in_band_branching(rep, band_ind) = ...
            sum((unique(final_grid_lineages_branching((band_ind+1):(band_ind+band_width),:))>0));

        num_spec_in_band_tube(rep, band_ind) = ...
            sum((unique(final_grid_lineages_tube((band_ind+1):(band_ind+band_width),:))>0));
         
         
        %loop over cells in band and compute average kappa_lin
        for cell_x = (band_ind+1):(band_ind+band_width)
            for cell_y = 1:grid_size_y


                %number of neighbour cells (fixed at 6)
                NUM_NEIGHBOURS = 6;
    

                %compute neighbours of this cell in the BRANCHING TREE
                [branch_neighbours_i, branch_neighbours_j] = get_neighbour_indices(cell_x,cell_y,...
                    size(final_grid_lineages_branching),generation_lengths,generation_circs);
    
                %get neighbour lineages
                this_cell_id = final_grid_lineages_branching(cell_x,cell_y);
                neighbour_ids = zeros(1,NUM_NEIGHBOURS);
                for neighbour_ind = 1:NUM_NEIGHBOURS
                    neighbour_ids(neighbour_ind) = final_grid_lineages_branching(...
                        branch_neighbours_i(neighbour_ind),...
                        branch_neighbours_j(neighbour_ind));
                end
    
                %compute proportion of same-lineage neighbours of this cell
                kappa_in_band_branching(rep, band_ind) = kappa_in_band_branching(rep, band_ind) + ...
                    sum(neighbour_ids==this_cell_id)/(NUM_NEIGHBOURS*band_width*grid_size_y);



                %compute neighbours of this cell in the TUBE
                [tube_neighbours_i, tube_neighbours_j] = get_neighbour_indices(cell_x,cell_y,...
                    size(final_grid_lineages_tube),sum(generation_lengths),generation_circs(1));
    
                %get neighbour lineages
                this_cell_id = final_grid_lineages_tube(cell_x,cell_y);
                neighbour_ids = zeros(1,NUM_NEIGHBOURS);
                for neighbour_ind = 1:NUM_NEIGHBOURS
                    neighbour_ids(neighbour_ind) = final_grid_lineages_tube(tube_neighbours_i(neighbour_ind),...
                        tube_neighbours_j(neighbour_ind));
                end
    
                %compute proportion of same-lineage neighbours of this cell
                kappa_in_band_tube(rep, band_ind) = kappa_in_band_tube(rep, band_ind) + ...
                    sum(neighbour_ids==this_cell_id)/(NUM_NEIGHBOURS*band_width*grid_size_y);
            end
        end
    end
end



%% plot

%make num species band plot
figure
hold on

%plot means
plot((1:num_bands)+1+band_width/2, mean(num_spec_in_band_branching,1), 'LineWidth', 1.5, 'Color', '#d95f02')
plot((1:num_bands)+1+band_width/2, mean(num_spec_in_band_tube,1), 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6])

%plot 10th to 90th percentiles
patch('XData', [(1:num_bands)+1+band_width/2, fliplr((1:num_bands)+1+band_width/2)],...
    'YData', [prctile(num_spec_in_band_branching,10,1), fliplr(prctile(num_spec_in_band_branching,90,1))], ...
    'FaceColor', '#d95f02', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch('XData', [(1:num_bands)+1+band_width/2, fliplr((1:num_bands)+1+band_width/2)],...
    'YData', [prctile(num_spec_in_band_tube,10,1), fliplr(prctile(num_spec_in_band_tube,90,1))], ...
    'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

%formatting
ylabel('Number of Lineages')
ylim([0,5])
yticks([0, 1, 2, 3, 4, 5])
xlabel('Depth (Cells)')
legend('Branching (Source Gen. 1)', 'Tube', 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 8)

%optionally save
if save_stuff
    SaveAsPngAndFig([], 'num_spec_plot', 9, 1.4, 9);
end



%make kappa_lin band plot
figure
hold on

%plot means
plot((1:num_bands)+1+band_width/2, mean(kappa_in_band_branching,1), 'LineWidth', 1.5, 'Color', '#d95f02')
plot((1:num_bands)+1+band_width/2, mean(kappa_in_band_tube,1), 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6])

%plot 10th to 90th percentiles
patch('XData', [(1:num_bands)+1+band_width/2, fliplr((1:num_bands)+1+band_width/2)],...
    'YData', [prctile(kappa_in_band_branching,10,1), fliplr(prctile(kappa_in_band_branching,90,1))], ...
    'FaceColor', '#d95f02', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch('XData', [(1:num_bands)+1+band_width/2, fliplr((1:num_bands)+1+band_width/2)],...
    'YData', [prctile(kappa_in_band_tube,10,1), fliplr(prctile(kappa_in_band_tube,90,1))], ...
    'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

%formatting
ylabel('Lineage Clusteredness $\kappa_{\mathrm{lin}}$')
xlabel('Depth (Cells)')
ylim([0.7, 1.1])
legend('Branching (Source Gen. 1)', 'Tube', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 8)

%optionally save
if save_stuff
    SaveAsPngAndFig([], 'kappa_lin_plot', 9, 1.4,9);
end



