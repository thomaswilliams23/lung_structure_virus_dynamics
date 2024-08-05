


%Makes Figure 5(a) from the manuscript


%% setup

%first set text intepreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath ../helper_functions
addpath plot_helpers


%data sources
open_edge_data = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_lineages/gen_1/';
branched_edge_data = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_lineages/gen_5/';


%tissue specs
num_reps = 100;
num_lineages = 4;


%define threshold for extinction
ext_thresh = 300;


%save stuff?
save_stuff = 0;


%% initialise

%initialise for lineage proportion counts
spec_num_edges = 0:0.1:1;
spec_num_centres = (spec_num_edges(1:end-1)+spec_num_edges(2:end))/2;

%initialise structure for number of cells for each species
prop_cells_given_species_gen_1 = zeros(1,num_lineages*num_reps);
prop_cells_given_species_gen_5 = zeros(1,num_lineages*num_reps);

%initialise probability of extinction for open edge case
prob_of_ext = 0;


%% loop

%loop iterations
for rep = 1:num_reps

    %load open edge data
    load(strcat(open_edge_data, '/sim_data_rep_',num2str(rep),'.mat'))
    open_edge_sim_data_this_rep = sim_data_this_rep;

    %load branched edge data
    load(strcat(branched_edge_data, '/sim_data_rep_',num2str(rep),'.mat'))
    branched_edge_sim_data_this_rep = sim_data_this_rep;


    %extract the final grid state and process it to replace every entry
    %with the lineage number
    mark_lineage_num = @(cell_id) get_infected_lineage_incl_dead(cell_id, num_lineages);
    
    open_edge_grid = open_edge_sim_data_this_rep.final_grid_state;
    open_edge_grid_species = arrayfun(mark_lineage_num, open_edge_grid);

    branched_edge_grid = branched_edge_sim_data_this_rep.final_grid_state;
    branched_edge_grid_species = arrayfun(mark_lineage_num, branched_edge_grid);


    %pull out lineage data
    num_ext_this_rep = 0;
    for lineage = 1:num_lineages

        %number of cells infected by this species
        prop_cells_given_species_gen_1((rep-1)*num_lineages + lineage) = ...
            sum(sum((open_edge_grid_species == lineage)))/numel(open_edge_grid_species);
        prop_cells_given_species_gen_5((rep-1)*num_lineages + lineage) = ...
            sum(sum((branched_edge_grid_species == lineage)))/numel(branched_edge_grid_species);

        %check for extinction
        if ~sum(sum((open_edge_grid_species((ext_thresh+1:end),:)==lineage)))
            num_ext_this_rep = num_ext_this_rep + 1;
        end
    end

    %update probability of extinction
    prob_of_ext = prob_of_ext + num_ext_this_rep/(num_reps*num_lineages);

end


%% plot histogram of lineage infected proportions


%first display probability of extinction
fprintf('\nProbability of extinction: %.3f\n', prob_of_ext);


%construct histogram
figure
hold on

[spec_num_counts_gen_1, spec_num_edges] = histcounts(prop_cells_given_species_gen_1, spec_num_edges);
[spec_num_counts_gen_5, spec_num_edges] = histcounts(prop_cells_given_species_gen_5, spec_num_edges);

b2 = bar(spec_num_centres, spec_num_counts_gen_5/sum(spec_num_counts_gen_5), 1, ...
    'FaceColor', '#66a61e', 'EdgeColor', 'none', 'FaceAlpha', 1);
b1 = bar(spec_num_centres, spec_num_counts_gen_1/sum(spec_num_counts_gen_1), 1, ...
    'FaceColor', '#d95f02', 'EdgeColor', '#d95f02', 'FaceAlpha', 0.2);

%formatting
xlim([0,1])
ylim([0,1])
xticks(spec_num_centres)
xticklabels({'0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', ...
    '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1'})
xlabel('Prop. Cell Sheet Infected by Species')
ylabel('Probability')
legend([b1, b2], {'Source Gen. 1', 'Source Gen. 5'}, 'interpreter', 'latex', 'FontSize', 8)


%optionally save
if save_stuff
    SaveAsPngAndFig(-1, 'lineage_histogram_branching_left_and_right_edge', 9, 1.4, 9);
end

