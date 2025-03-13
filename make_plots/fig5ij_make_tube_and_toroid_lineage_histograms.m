


%Makes Figures 5(i) and (j) from the manuscript


%% setup

%first set text intepreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath plot_helpers


%data sources
tube_data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_tube_lineages_full_data/tube_circ_';
toroid_data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_toroid_lineages_full_data/tube_circ_';


%loop specs
all_circs = [4,8,16,32,64,128];
circ_labels = {'4 cells', '8 cells', '16 cells', '32 cells', '64 cells', '128 cells'};


%tissue specs
num_reps = 100;
num_lineages = 4;


%define threshold for extinction
ext_thresh = 300;


%colours
tube_colours = winter(length(all_circs)+1);


%save stuff?
save_stuff = 1;



%% initialise

%initialise for lineages proportion counts
spec_prop_edges = 0:0.1:1;
spec_prop_centres = (spec_prop_edges(1:end-1)+spec_prop_edges(2:end))/2;
all_lin_props_tube = zeros(length(all_circs), num_lineages*num_reps);
all_lin_props_toroid = zeros(length(all_circs), num_lineages*num_reps);

%initialise for species progression depth
all_inv_depths = zeros(length(all_circs), num_lineages*num_reps);

%initialise for extinction
prob_of_ext = zeros(1,length(all_circs));



%% loop

%loop circumferences
for circ_ind = 1:length(all_circs)

    circ=all_circs(circ_ind);


    %initialise structure for all positions of lineage dieout
    last_col_this_circ = zeros(1,num_lineages*num_reps);

    %initialise structure for number of cells for each lineage
    prop_cells_given_species_toroid = zeros(1,num_lineages*num_reps);
    prop_cells_given_species_tube = zeros(1,num_lineages*num_reps);

    %initialise count of species which go extinct
    num_ext = 0;


    %loop iterations
    for rep = 1:num_reps


        %load toroid data
        load(strcat(toroid_data_source, num2str(circ), '/sim_data_',num2str(rep),'.mat'))
        toroid_sim_data_this_rep = sim_data_this_rep;

        %load tube data
        load(strcat(tube_data_source, num2str(circ), '/sim_data_',num2str(rep),'.mat'))
        tube_sim_data_this_rep = sim_data_this_rep;


        %extract the final grid states and process to replace every entry
        %with the lineage number
        mark_lineage_num = @(cell_id) get_infected_lineage_incl_dead(cell_id, num_lineages);

        toroid_cell_grid =toroid_sim_data_this_rep.final_grid_state;
        toroid_cell_grid_species = arrayfun(mark_lineage_num,toroid_cell_grid);
        
        tube_cell_grid = tube_sim_data_this_rep.final_grid_state;
        tube_cell_grid_species = arrayfun(mark_lineage_num, tube_cell_grid);


        %pull out lineage data
        ext_this_rep = 0;
        for lineage = 1:num_lineages

            %count up how many of this lineage at each column and find the
            %last one
            num_lin_at_col = sum((tube_cell_grid_species == lineage),2);
            last_col = find(num_lin_at_col,1,'last');

            %save out and increment extinction instances
            last_col_this_circ((rep-1)*num_lineages + lineage) = last_col;
            num_ext = num_ext + (last_col<ext_thresh);

            %proportion of cells infected by this lineage
            prop_cells_given_species_tube((rep-1)*num_lineages + lineage) = ...
                sum(sum((tube_cell_grid_species == lineage)))/numel(tube_cell_grid_species);
            prop_cells_given_species_toroid((rep-1)*num_lineages + lineage) = ...
                sum(sum((toroid_cell_grid_species == lineage)))/numel(toroid_cell_grid_species);
        end

    end

    %save out all histogram counts and invasion depths
    all_lin_props_tube((all_circs==circ), :) = prop_cells_given_species_tube;
    all_lin_props_toroid((all_circs==circ),:) = prop_cells_given_species_toroid;
    all_inv_depths((all_circs==circ),:) = last_col_this_circ;


    %plot histogram of lineage infected proportions
    figure
    hold on

    [spec_num_counts_tube, spec_prop_edges] = histcounts(prop_cells_given_species_tube, spec_prop_edges);
    [spec_num_counts_toroid, spec_prop_edges] = histcounts(prop_cells_given_species_toroid, spec_prop_edges);

    bar(spec_prop_centres, spec_num_counts_toroid/sum(spec_num_counts_toroid), 1, ...
        'FaceColor', [0.7,0.7,0.7], 'EdgeColor', 'none');
    bar(spec_prop_centres, spec_num_counts_tube/sum(spec_num_counts_tube), 1, ...
        'FaceColor',tube_colours(circ_ind,:), 'FaceAlpha', 0.8);


    %formatting
    xlim([0,1])
    ylim([0,1])

    xticks(spec_prop_centres)
    xticklabels({'0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', ...
        '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1'})

    xlabel('Prop. Cell Sheet Infected by Lineage')
    ylabel('Probability')

    legend({'Toroid', 'Tube'}, 'interpreter', 'latex', 'FontSize', 8)


    %save plot
    if save_stuff
        SaveAsPngAndFig(-1, strcat('lineage_histogram_circ_',num2str(circ)), 8, 1.4, 9);
    end


    %update probability of extinction
    prob_of_ext((all_circs==circ)) = num_ext/(num_reps*num_lineages);

end

%display probability of extinction
fprintf('\nProbabilities of extinction:\n Circs:\n')
disp(all_circs)
fprintf('\nProb ext:\n')
disp(prob_of_ext)