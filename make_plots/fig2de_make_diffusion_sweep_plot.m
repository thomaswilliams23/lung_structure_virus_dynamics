


%Makes Figure 2(d) and (e) from the manuscript



%% data specs and initialisation

%first set text interpreter
set(0, 'DefaultTextInterpreter', 'latex')

%and add path
addpath plot_helpers




%data source
source_folder = '../SIM_DATA_from_ms/SWEEP_OUTPUT_diffusion_and_aspect_ratio_sweep/';
template_1 = '/tube_circ_';
template_2 = '_diff_';

%tube circumferences
circ_folders = {'4', '8', '16', '32', '64'};
circ_vals = 2*2.^(1:5);

%diffusion coefficients
diff_folders = {'0', '1', '10', '100', '1000', 'Inf'};
diff_legend = {'0', '1', '10', '100', '1000', '$\infty$'};
diff_vals_all = [0, 1, 10, 100, 1000, Inf];
diff_vals_internal = diff_vals_all(2:end-1)';

%tube lengths and aspect ratios
tot_cells = 4096;
tubes_lengths = tot_cells./circ_vals;
AR_vals = circ_vals./tubes_lengths;

%initialise s_D
t_0_over_t_inf = zeros(1,length(circ_vals));
t_100_over_t_inf = zeros(1,length(circ_vals));

%number iterations per AR-diffusion combination
num_reps = 10;

%time specs
dt = 0.1;
def_max_num_time_pts = 60000;

%colours
diff_plot_colours = min(1.3*summer(length(diff_folders)+3),1);

%save output?
save_stuff = 1;



%% main loop
    
%loop over AR
for circ_ind = 1:length(circ_folders)

    %initialise
    inv_times_this_circ = zeros(length(diff_folders), num_reps);
    max_num_time_pts = 0;

    %set up plot
    diff_plot_handles = zeros(length(diff_folders),1);
    figure(circ_ind)
    hold on


    %loop over diffusion
    for diff_ind = 1:length(diff_folders)

        %load data
        load(strcat(source_folder, template_1, ...
            circ_folders{circ_ind}, template_2, ...
            diff_folders{diff_ind}, '/all_sim_data_this_fldr'));

        %initialise an array for time to (95%) invasion
        inv_times_this_diff = zeros(1,num_reps);


        %loop over reps and find invasion times
        for rep = 1:num_reps

            DIEOUT_FLAG = -1;
            invaded_thresh = 0.95;

            if all_sim_data_this_fldr{rep}.prop_infected(end)>invaded_thresh
                final_ind = find(all_sim_data_this_fldr{rep}.prop_infected>invaded_thresh, 1, 'first');
                inv_times_this_diff(rep) = dt*final_ind;
            else
                inv_times_this_diff(rep) = DIEOUT_FLAG;
            end


            %add to the plot
            diff_plot_handles(diff_ind) = plot(dt*(0:all_sim_data_this_fldr{rep}.num_output_points-1),...
                    all_sim_data_this_fldr{rep}.prop_infected, ...
                    'color', diff_plot_colours(diff_ind,:), 'LineWidth', 0.5);

            %account for possibility of early termination of simulation
            if all_sim_data_this_fldr{rep}.num_output_points < def_max_num_time_pts
                max_num_time_pts = max(max_num_time_pts, all_sim_data_this_fldr{rep}.num_output_points);
            end

        end

        %save out
        inv_times_this_circ(diff_ind,:) = inv_times_this_diff;

    end


    %add annotations
    plot([0, 1.1*max_num_time_pts*dt], [1,1], 'k--')
    text(0.8*max_num_time_pts*dt, 1.05, 'Tissue destroyed', ...
        'Interpreter', 'latex', 'FontSize',8)

    plot([0, 1.1*max_num_time_pts*dt], 0.95*[1,1], 'r--', 'LineWidth', 1.2)
    text(0.66*max_num_time_pts*dt, 0.86, '\bf{95\% tissue destroyed}', ...
        'Interpreter','latex', 'Color', 'r', 'FontSize',8)


    %add dropdowns
    for diff_ind = 1:size(inv_times_this_circ,1)
        time_to_inv_this_diff = inv_times_this_circ(diff_ind,:);
        mean_t_inv = mean(time_to_inv_this_diff(time_to_inv_this_diff>0));
        if diff_ind == 1 || diff_ind == size(inv_times_this_circ,1)
            plot([mean_t_inv, mean_t_inv], [0, 0.95], 'r:', 'LineWidth', 1.1)
        else
            plot([mean_t_inv, mean_t_inv], [0, 0.95], 'r:', 'LineWidth', 0.5)
        end
    end


    %formatting
    legend(diff_plot_handles, diff_legend, 'Location', 'EastOutside', 'Interpreter', 'latex');
    xlim([0, 1.1*max_num_time_pts*dt])
    ylim([0, 1.2])
    xlabel('Time (h)')
    ylabel('Cumul. Proportion Infected')


    %save
    if save_stuff
        SaveAsPngAndFig([], strcat('cumul_time_series_circ_', circ_folders{circ_ind}), 12, 2, 9);
    end




    %compute t_0/t_inf and t_100/t_inf
    t_0_over_t_inf(circ_ind) = inv_times_this_circ(1)/inv_times_this_circ(end);
    t_100_over_t_inf(circ_ind) = inv_times_this_circ(4)/inv_times_this_circ(end);

end



% make a plot of t_0/t_inf and t_100/t_inf vs AR
figure

diff_sens_plt = semilogx(AR_vals, t_0_over_t_inf, 'color', [0.8500 0.3250 0.0980]);
set(diff_sens_plt, 'linewidth', 1.5)
set(diff_sens_plt, 'marker', 'o')

hold on

diff_sens_plt_2 = semilogx(AR_vals, t_100_over_t_inf, 'color', [0.4940 0.1840 0.5560]);
set(diff_sens_plt_2, 'linewidth', 1.5)
set(diff_sens_plt_2, 'marker', 'o')

legend([diff_sens_plt, diff_sens_plt_2], {'$t_{0.95}^{0}/t_{0.95}^{\infty}$',...
    '$t_{0.95}^{100}/t_{0.95}^{\infty}$'}, 'Interpreter', 'latex', 'FontSize', 9, ...
    'Location', 'northeast')

xlim([10^-2.6, 10^0.2])
xlabel('Aspect Ratio')
ylabel('Fold Increase in Infection Time')

box off



%save
if save_stuff
    SaveAsPngAndFig([], 'rate_ratio_plot', 8, 4/3, 9);
    ax = gca;
    ax.Title.FontSize = 11;
    exportgraphics(gcf, 'rate_ratio_plot.png', 'Resolution',300)
end


