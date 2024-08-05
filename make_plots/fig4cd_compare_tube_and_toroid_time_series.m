


%Makes Figure 4(c) and (d) from the manuscript



%% setup

%first set text intepreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath plot_helpers


%data sources
tube_data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_tube_lineages_for_time_series/tube_circ_';
toroid_data_source = '../SIM_DATA_from_ms/SWEEP_OUTPUT_toroid_lineages_for_time_series/tube_circ_';


%loop specs
all_circs = [4,8,16,32,64,128];
num_reps_to_plot = 10;


%time
dt = 0.5;
num_output_points = 400;
time_vec = dt*(0:(num_output_points-1));


%colours (for tubes)
plot_colours = winter(length(all_circs)+1);

%save stuff?
save_stuff = 0;


%% loop folders
for circ = all_circs


    %initialise
    figure 
    hold on

    mean_toroid_time_half_inv = 0;
    mean_tube_time_half_inv = 0;


    %% toroidal data
    num_non_dieout = 0;

    %loop iterations
    for rep = 1:num_reps_to_plot

        %load data
        load(strcat(toroid_data_source, num2str(circ), '/sim_data_', num2str(rep)));

        %add to plot
        p1 = plot(time_vec, sim_data_this_rep.net_I, ...
            'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.2);

        %compute time to half infection
        toroid_half_inv_ind = time_to_half_inv(sim_data_this_rep.net_I);

        %check for dieout case and update mean time to half infection
        TIME_HALF_INV_DIEOUT_CUTOFF = 10;
        if toroid_half_inv_ind*dt>TIME_HALF_INV_DIEOUT_CUTOFF
            num_non_dieout = num_non_dieout + 1;
            mean_toroid_time_half_inv = mean_toroid_time_half_inv + toroid_half_inv_ind*dt;
        end

    end

    %compute mean over non-dieout cases
    mean_toroid_time_half_inv = mean_toroid_time_half_inv/num_non_dieout;



    %% tube data
    num_non_dieout = 0;

    %loop iterations
    for rep = 1:num_reps_to_plot

        %load data
        load(strcat(tube_data_source, num2str(circ), '/sim_data_', num2str(rep)));

        %add to plot
        p2 = plot(time_vec, sim_data_this_rep.net_I, ...
            'Color', plot_colours((all_circs==circ),:), 'LineWidth', 0.2);

        %compute time to half infection
        tube_half_inv_ind = time_to_half_inv(sim_data_this_rep.net_I);

        %check for dieout case and update mean time to half infection
        TIME_HALF_INV_DIEOUT_CUTOFF = 10;
        if tube_half_inv_ind*dt>TIME_HALF_INV_DIEOUT_CUTOFF
            num_non_dieout = num_non_dieout + 1;
            mean_tube_time_half_inv = mean_tube_time_half_inv + tube_half_inv_ind*dt;
        end

    end

    %compute mean over non-dieout cases
    mean_tube_time_half_inv = mean_tube_time_half_inv/num_non_dieout;



    %% time to half infection  

    %work out height on curves at time to half infection
    max_toroid_half_inv_height = 0;
    max_tube_half_inv_height = 0;

    %toroid
    for rep = 1:num_reps_to_plot
        load(strcat(toroid_data_source, num2str(circ), '/sim_data_', num2str(rep)));
        max_toroid_half_inv_height = max(max_toroid_half_inv_height, sim_data_this_rep.net_I(round(mean_toroid_time_half_inv/dt)));
    end

    %tube
    for rep = 1:num_reps_to_plot
        load(strcat(tube_data_source, num2str(circ), '/sim_data_', num2str(rep)));
        max_tube_half_inv_height = max(max_tube_half_inv_height, sim_data_this_rep.net_I(round(mean_tube_time_half_inv/dt)));
    end


    %plot time to half infection
    plot(mean_toroid_time_half_inv*[1,1], [0, max_toroid_half_inv_height], ...
        'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.8, 'LineStyle', '--')
    plot(mean_tube_time_half_inv*[1,1], [0, max_tube_half_inv_height], ...
        'Color', plot_colours((all_circs==circ),:), 'LineWidth', 0.8, 'LineStyle', '--')


    %annotations
    text(80, 0.3, sprintf('$t_{0.5} = %.1f$ h', mean_toroid_time_half_inv), ...
        'FontSize', 10, 'Color', [0.3, 0.3, 0.3], 'Interpreter', 'latex')
    text(120, 0.15, sprintf('$t_{0.5} = %.1f$ h', mean_tube_time_half_inv), ...
        'FontSize', 10, 'Color', plot_colours((all_circs==circ),:), 'Interpreter','latex')


    %formatting
    xlim([0,200])
    ylim([0,0.5])
    xlabel('Time (h)')
    ylabel('Prop. Infected')
    legend([p1, p2], {'Toroid', 'Tube'}, 'Interpreter', 'latex', 'FontSize', 8)


    %save?
    if save_stuff
        mkdir('time_series_tubes_vs_toroids')
        SaveAsPngAndFig([], strcat('time_series_tubes_vs_toroids/circ_',num2str(circ)),8, 1.6, 9);
    end


end

