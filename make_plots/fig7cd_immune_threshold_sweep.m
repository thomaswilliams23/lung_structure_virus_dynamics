


%Makes Figure 7(c) and (d) from the manuscript



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
num_reps = 100;

%time
max_t = 300;
dt = 0.1;
num_time_points = round(max_t/dt);

%immune parameters
F_thresh_sweep = 0.1:0.05:0.5;
F_thresh_labels = {};
for F_thresh_ind = 1:length(F_thresh_sweep)
    F_thresh_labels{F_thresh_ind} = num2str(F_thresh_sweep(F_thresh_ind));
end
t_act = 30;
t_act_timesteps = round(t_act/dt);

%legend
legend_names = {'Mean $F^{\infty}$', '$F^{\infty}$ of mean trajectory'};

%colours
open_edge_colour = [217, 95, 2]/256; 
branching_edge_colour = [102, 166, 30]/256;

%save?
save_stuff = 1;





%% plot
for start_on_left_edge = [0,1]


    %initialise
    av_time_series = zeros(1,num_time_points);
    num_non_dieout = 0;
    F_inf_all = zeros(num_reps, length(F_thresh_sweep));
    F_inf_av_time_series = zeros(1, length(F_thresh_sweep));


    %loop reps
    for rep = 1:num_reps

        %load in data
        if start_on_left_edge
            load(strcat(open_edge_data, '/sim_data_', num2str(rep), '.mat'));
        else
            load(strcat(branching_edge_data, '/sim_data_', num2str(rep), '.mat'));
        end
        this_sim = sim_data_this_rep;


        %filter out dieout cases
        if this_sim.prop_infected(end)<max(F_thresh_sweep)
            F_inf_all(rep, :) = NaN;
            continue
        end

        %update averages
        num_non_dieout = num_non_dieout + 1;
        av_time_series = av_time_series + this_sim.prop_infected;


        %loop immune thresholds
        for F_thresh_ind = 1:length(F_thresh_sweep)

            %set immune threshold
            F_thresh = F_thresh_sweep(F_thresh_ind);

            %find index of threshold and immune activation
            ind_at_thresh = find(this_sim.prop_infected>F_thresh, 1, 'first');
            ind_at_act = ind_at_thresh+t_act_timesteps;

            %case where activation time exceeds final time (tissue
            %destroyed)
            if ind_at_act>num_time_points
                F_inf_all(rep, F_thresh_ind) = 1;
            %otherwise test damage at activation time
            else
                F_inf_all(rep, F_thresh_ind) = this_sim.prop_infected(ind_at_act);
            end


        end

    end
    
    %trim nans (dieouts)
    non_nan_cols = find(~isnan(F_inf_all(:,1)));
    F_inf_all = F_inf_all(non_nan_cols, :);

    %update av time series
    av_time_series = av_time_series/num_non_dieout;


    %now work out F_inf on averaged time series
    for F_thresh_ind = 1:length(F_thresh_sweep)

        %set immune threshold
        F_thresh = F_thresh_sweep(F_thresh_ind);

        %find index of threshold and immune activation
        ind_at_thresh = find(av_time_series>F_thresh, 1, 'first');
        ind_at_act = ind_at_thresh+t_act_timesteps;

        %case where activation time exceeds final time (tissue
        %destroyed)
        if ind_at_act>num_time_points
            F_inf_av_time_series(F_thresh_ind) = 1;
        %otherwise test damage at activation time
        else
            F_inf_av_time_series(F_thresh_ind) = av_time_series(ind_at_act);
        end


    end



    %% plot

    figure
    hold on

    %retrieve colour
    if start_on_left_edge
        colour_this_rep = open_edge_colour;
    else
        colour_this_rep = branching_edge_colour;
    end

    %violins
    [~,~,~,~,~,av_of_f,kern_dens,kern_x]=violin(F_inf_all, 'facecolor', colour_this_rep, ...
        'edgecolor', colour_this_rep, 'facealpha', 0.5, 'medc', '');
    legend off
    
    %F_inf on mean time series
    plot(1:length(F_thresh_sweep), F_inf_av_time_series, 'r-', 'LineWidth', 0.5);

    %plot horizontals on violins
    for vln_ind = 1:length(F_thresh_sweep)
        kern_width = interp1(kern_x(:,vln_ind),kern_dens(:,vln_ind), F_inf_av_time_series(vln_ind));
        f_of_av= plot([vln_ind-kern_width, vln_ind+kern_width], F_inf_av_time_series(vln_ind)*[1,1], 'r-', 'LineWidth', 1);
    end


    %% formatting

    if start_on_left_edge
        title_text = 'Source Open Edge';
        fname = 'open_edge';
    else
        title_text = 'Source Branched Edge';
        legend([av_of_f, f_of_av], legend_names, 'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 8);
        fname = 'branching_edge';
    end

    text(0.7, 0.85, title_text, 'FontSize', 11, 'Interpreter','latex')

    xticks(1:length(F_thresh_sweep))
    xticklabels(F_thresh_labels)

    ylim([0, 0.9])

    xlabel('$F^{\mathrm{thresh}}$')
    ylabel('$F^{\infty}$')

    box off
    

    %optionally save
    if save_stuff
        SaveAsPngAndFig([], strcat('immune_threshold_sweep_', fname), 9, 4/3, 9);
    end

end




