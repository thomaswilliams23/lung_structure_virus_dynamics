


%Makes Figure 3(c) from the manuscript



%% setup

%first set text interpreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath plot_helpers


%data source
base_folder = '../SIM_DATA_from_ms/SWEEP_OUTPUT_generation_seed_sweep/';
folder_template = 'x_seed_position_';


%loop specs
num_gens = 5;
gen_lims = cumsum([100, 100, 100, 100, 100]);
folder_indices_to_plot = [1, 101, 201, 301, 401];
num_reps = 10;


%time
dt = 0.1;
max_time = 600;
max_time_to_plot = 300;
max_num_points = round(max_time/dt);


%colours
plot_colours = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a",  "#66a61e"];


%save?
save_stuff = 1;


%% initialise
plt_handles_net_I = 0*gen_lims;

%means
mean_net_I = zeros(num_gens,max_num_points);
max_num_valid_points = max_num_points*ones(1,num_gens);


%figure
figure
hold on


%% loop over folders
for source_gen = length(gen_lims):-1:1

    %find folder
    x_seed_position = folder_indices_to_plot(source_gen);
    source_folder = strcat(base_folder, folder_template, num2str(x_seed_position));


    %load data
    load(strcat(source_folder, '/sim_data_this_seed_position'));


    %loop iterations
    for iter = 1:num_reps

        %extract simulation
        this_sim = sim_data_this_seed_position{iter};

        %plot
        plt_handles_net_I(source_gen) = patchline(dt*(0:(this_sim.num_output_points-1)), ...
            this_sim.net_I,... %/this_sim.prop_infected(end)
            'edgecolor', plot_colours(source_gen), 'LineWidth',0.3, 'EdgeAlpha', 0.3);


        %tally mean
        mean_net_I(source_gen, 1:this_sim.num_output_points) = ...
            mean_net_I(source_gen, 1:this_sim.num_output_points) + this_sim.net_I/num_reps;


        %this keeps track of the maximum number of output points before
        %simulation terminated
        max_num_valid_points(source_gen) = min(max_num_valid_points(source_gen), this_sim.num_output_points);
    end
end


%% plot averages
for source_gen = 1:num_gens

    plt_handles_net_I(source_gen) = plot(dt*(0:max_num_valid_points(source_gen)-1), ...
        mean_net_I(source_gen, 1:max_num_valid_points(source_gen)), ...
        'color', plot_colours(source_gen), 'LineWidth', 1.5);

end


%% formatting
legend(plt_handles_net_I, {'Gen. 1', 'Gen. 2', 'Gen. 3', 'Gen. 4', 'Gen. 5'}, ...
    'Location', 'NorthEast','Interpreter', 'latex')
ylabel('Instantaneous Infected Proportion')
ylim([0,0.2])
xlim([0,max_time_to_plot])
xlabel('Time (h)')

%optionally save
if save_stuff
    SaveAsPngAndFig([], 'branching_tube_net_dynamics', 9, 1.4, 9);
end