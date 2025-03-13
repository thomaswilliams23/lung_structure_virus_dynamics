


%Makes Figure 8(b) from the manuscript



%% setup

%first set text intepreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath ../helper_functions
addpath plot_helpers


%data sources
open_edge_data = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_top_and_bottom_seed/open_edge';
branching_edge_data = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_top_and_bottom_seed/branching_edge';


%number time series to plot
num_reps_to_plot = 10;

%time
max_t = 300;
dt = 0.1;
num_time_points = round(max_t/dt);
time_vec = dt*(0:(num_time_points-1));

%legend
legend_names = {'Open Edge', 'Branched Edge'};

%colours
open_edge_colour = [217, 95, 2]/256; 
branching_edge_colour = [102, 166, 30]/256;

%save?
save_stuff = 1;



%% plot

%set up figure
figure
hold on

%initialise
av_time_series_open_edge = zeros(1,num_time_points);
num_non_dieout_open_edge = 0;

av_time_series_branching_edge = zeros(1,num_time_points);
num_non_dieout_branching_edge = 0;


%loop reps
for rep = 1:num_reps_to_plot

    %load in data
    load(strcat(open_edge_data, '/sim_data_', num2str(rep), '.mat'));
    this_sim_open_edge = sim_data_this_rep;

    load(strcat(branching_edge_data, '/sim_data_', num2str(rep), '.mat'));
    this_sim_branching_edge = sim_data_this_rep;

    %plot
    patchline(time_vec, this_sim_open_edge.prop_infected, ...
        'edgecolor', open_edge_colour, 'LineWidth', 0.3, 'EdgeAlpha', 0.4);

    patchline(time_vec, this_sim_branching_edge.prop_infected, ...
        'edgecolor', branching_edge_colour, 'LineWidth', 0.3, 'EdgeAlpha', 0.4);


    %add to averages if non-dieout
    DIEOUT_CUTOFF = 0.9;
    if this_sim_open_edge.prop_infected(end)>DIEOUT_CUTOFF
        num_non_dieout_open_edge = num_non_dieout_open_edge + 1;
        av_time_series_open_edge = av_time_series_open_edge + this_sim_open_edge.prop_infected;
    end
    if this_sim_branching_edge.prop_infected(end)>DIEOUT_CUTOFF
        num_non_dieout_branching_edge = num_non_dieout_branching_edge + 1;
        av_time_series_branching_edge = av_time_series_branching_edge + this_sim_branching_edge.prop_infected;
    end

end

%update averages
av_time_series_open_edge = av_time_series_open_edge/num_non_dieout_open_edge;
av_time_series_branching_edge = av_time_series_branching_edge/num_non_dieout_branching_edge;


%add to plot
p1 = plot(time_vec, av_time_series_open_edge,...
    'Color', open_edge_colour, 'LineWidth', 1.5);
p2 = plot(time_vec, av_time_series_branching_edge,...
    'Color', branching_edge_colour, 'LineWidth', 1.5);



%also plot tissue destroyed line
plot([0, max_t], [1,1], 'k--')
text(max_t-100, 1+0.05, 'Tissue Destroyed', 'Interpreter', 'latex', 'FontSize', 8)



%% formatting

xlim([0, max_t]);
ylim([0, 1.2]);
xlabel('Time (h)')
ylabel('Cumulative Infected Proportion')
legend([p1, p2], legend_names, ...
    'Interpreter', 'latex', 'FontSize', 8, 'Location', 'SouthEast')


%optionally save
if save_stuff
    SaveAsPngAndFig([], 'time_series_open_and_branching_edge', 9, 4/3, 9);
end
