


%Makes Figure 2(d) from the manuscript



%% setup

%first set text interpreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath plot_helpers


%data source
base_folder = '../SIM_DATA_from_ms/SWEEP_OUTPUT_generation_seed_sweep/';
folder_template = 'x_seed_position_';
folder_indices_to_plot = [1, 101, 201, 301, 401];

%loop specs
num_gens = 5;
gen_lims = cumsum([100,100,100,100,100]);
num_reps = 10;


%time
dt = 0.1;
max_time = 600;
max_num_points = round(max_time/dt);
max_plt_time = 300;
time_vec = dt*(0:(max_num_points-1));


%colours
plot_colours = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a",  "#66a61e"];


%save?
save_stuff = 1;



%% loop over seed cells positions
for fldr_to_plot = folder_indices_to_plot

    %load data
    source_folder = strcat(base_folder, folder_template, num2str(fldr_to_plot));
    load(strcat(source_folder, '/sim_data_this_seed_position.mat'))


    %% initialise

    %set up figure
    plt_handles_net_I = 0*gen_lims;
    figure

    %means
    mean_net_I = zeros(num_gens+1,max_num_points);
    max_num_valid_points = max_num_points;

    mean_peak_time = 0;
    mean_time_half_inv = 0;


    %% loop over iterations
    num_non_dieout = 0;
    for iter = 1:num_reps

        %load data
        sim_data_this_rep = sim_data_this_seed_position{iter};


        %update mean peak time and median - don't count dieout iterations
        [peak_height,peak_ind] = max(sim_data_this_rep.net_I);
        DIEOUT_CUTOFF = 0.05;

        this_inst_non_dieout = (peak_height>DIEOUT_CUTOFF);
        if this_inst_non_dieout
            num_non_dieout = num_non_dieout + 1;
            mean_peak_time = mean_peak_time + peak_ind*dt;
            mean_time_half_inv = mean_time_half_inv + time_to_half_inv(sim_data_this_rep.net_I)*dt;
        end



        %now plot
        for gen = 1:num_gens
            patchline(time_vec, sim_data_this_rep.section_sums(gen,:), ...
                'edgecolor', plot_colours(gen), 'LineWidth',0.3, 'EdgeAlpha', 0.3);
            hold on
        end
        patchline(time_vec, sim_data_this_rep.net_I, 'edgecolor', 'k', 'LineWidth',0.3, 'EdgeAlpha', 0.3);


        %update mean
        if this_inst_non_dieout
            for gen = 1:num_gens
                mean_net_I(gen, :) = mean_net_I(gen, :) + sim_data_this_rep.section_sums(gen,:);
            end
            mean_net_I(num_gens+1, :) = mean_net_I(num_gens+1, :) + sim_data_this_rep.net_I;
        end

    end

    %compute means over non-dieout instances
    mean_peak_time = mean_peak_time/num_non_dieout;
    mean_net_I = mean_net_I/num_non_dieout;
    mean_time_half_inv = mean_time_half_inv/num_non_dieout;


    %plot means
    for gen = 1:num_gens
        plt_handles_net_I(gen) = plot(time_vec, mean_net_I(gen,:),...
            'color', plot_colours(gen), 'LineWidth', 1.5);
    end
    plt_handles_net_I(num_gens+1) = plot(time_vec, mean_net_I(num_gens+1,:),...
        'color', 'k', 'LineWidth', 1.5);


    %plot peak
    peak_height = mean_net_I(num_gens+1,round(mean_peak_time/dt));
    plot(mean_peak_time*[1,1], [0, peak_height],...
        'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
    text(mean_peak_time-40, peak_height + 0.03, strcat(['$t_{\mathrm{peak}}$ = ', sprintf('%.1f h', mean_peak_time)]), ...
        'FontSize', 9, 'Color', 'k', 'Interpreter', 'latex')

    %plot time to half infected
    time_to_half_inv_height = mean_net_I(num_gens+1,round(mean_time_half_inv/dt));
    plot(mean_time_half_inv*[1,1], [0, time_to_half_inv_height],...
        'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    text(mean_time_half_inv-40, time_to_half_inv_height + 0.03, sprintf('$t_{0.5} = %.1f$ h', mean_time_half_inv), ...
        'FontSize', 9, 'Color', 'r', 'Interpreter', 'latex')


    %find source generation
    source_gen = find(folder_indices_to_plot==fldr_to_plot, 1, 'first');


    %formatting
    %legend(plt_handles_net_I, {'Gen. 1', 'Gen. 2', 'Gen. 3', 'Gen. 4', 'Gen. 5', 'Overall'},...
    %    'Location', 'EastOutside', 'Interpreter', 'latex')
    ylabel('Inst. Inf. Prop.')
    xlabel('Time (h)')
    xlim([0,max_plt_time])
    ylim([0, 0.2])
    title(sprintf('Source Gen. %d', source_gen))


    
    %optionally save
    if save_stuff
        SaveAsPngAndFig([], strcat('section_dynamics/section_dynamics_source_gen_', ...
            num2str(source_gen)), 6, 1.4, 9);
    end


end