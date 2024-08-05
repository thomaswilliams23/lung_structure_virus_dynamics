


%Makes Figure 3 from the manuscript



%% setup

%first set text interpreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath plot_helpers


%data source
base_folder = '../SIM_DATA_from_ms/SWEEP_OUTPUT_fine_res_seed_sweep/';
folder_template = 'x_seed_position_';
folder_inds = 1:25:500;


%loop specs
num_gens = 5;
gen_lims = cumsum([100, 100, 100, 100, 100]);
num_reps = 10;


%time
max_t = 240;
dt = 0.1;
plot_max = max_t;


%colours
gen_colours = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a",  "#66a61e"];
tpeak_colour = [0,0,0];
t_05_colour = [1,0,0];
dieout_colour = [1,0,0];


%font and format
fig_font_size = 9;
fig_width = 18;
fig_ratio = 3;


%save?
save_stuff = 1;


%% set up figure

figure
hold on

%loop generations
for i=1:length(gen_lims)
    plot(gen_lims(i)*[1,1], [0,plot_max], 'k--')

    if i==1

        %shade by generation colour
        patch('XData', [0, gen_lims(i), gen_lims(i), 0],...
            'YData', [0, 0, plot_max, plot_max], ...
            'FaceColor', gen_colours(i), 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

        %label generation
        text(0.34*gen_lims(i), 0.9*plot_max, sprintf('Gen. %d',i), 'FontSize', fig_font_size)
    else

        %shade by generation colour
        patch('XData', [gen_lims(i-1), gen_lims(i), gen_lims(i), gen_lims(i-1)],...
            'YData', [0, 0, plot_max, plot_max], ...
            'FaceColor', gen_colours(i), 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

        %label generation
        text(gen_lims(i-1)+0.34*(gen_lims(i)-gen_lims(i-1)), 0.9*plot_max, sprintf('Gen. %d',i), 'FontSize', fig_font_size)
    end

end


%plot a line for dieout cases
plot([0, gen_lims(end)], [0,0], 'r:', 'LineWidth', 1.2)
text(10, 0.05*plot_max, 'Dieout', 'Color', 'r', 'FontSize', round(0.9*fig_font_size))



%% loop seed positions

%initialise array for means at seed position
peak_time_means = zeros(1, length(folder_inds));
time_half_inv_means = zeros(1, length(folder_inds));

%loop
for start_pos_ind = 1:length(folder_inds)

    %load data
    target_folder = strcat(base_folder, folder_template, num2str(folder_inds(start_pos_ind)));
    load(strcat(target_folder, '/sim_data_this_seed_position'));

    %initialise arrays for this seed position
    time_to_peak_all = zeros(num_reps, 1);
    time_half_inv_all = zeros(num_reps, 1);
    sim_not_dieout = zeros(num_reps, 1);


    %compute time to reach threshold
    for rep = 1:num_reps

        this_sim = sim_data_this_seed_position{rep};

        %check for dieout
        not_dieout = (max(this_sim.net_I)>10*this_sim.net_I(1));
        sim_not_dieout(rep) = not_dieout;

        %record peak time
        if not_dieout
            [~,peak_ind] = max(this_sim.net_I);
            time_to_peak_all(rep) = peak_ind*dt;
            time_half_inv_all(rep) = time_to_half_inv(this_sim.net_I)*dt;
        end

    end

    %determine which simulations are not dieout cases
    which_sims_valid = find(sim_not_dieout);
    which_sims_dieout = find(~sim_not_dieout);


    %add to the plot
    if sum(which_sims_valid)

        %scatter t_peak and time to half infection
        scatter(folder_inds(start_pos_ind)*ones(size(time_to_peak_all(which_sims_valid))),...
            time_to_peak_all(which_sims_valid), ...
            50, tpeak_colour, '.', 'LineWidth', 2);
        scatter(folder_inds(start_pos_ind)*ones(size(time_half_inv_all(which_sims_valid))),...
            time_half_inv_all(which_sims_valid), ...
            50, t_05_colour, '.', 'LineWidth', 2);

        %compute means
        peak_time_means(start_pos_ind) = mean(time_to_peak_all(which_sims_valid));
        time_half_inv_means(start_pos_ind) = mean(time_half_inv_all(which_sims_valid));


    %error case (all dieout)
    else
        peak_time_means(start_pos_ind) = NaN;
        time_half_inv_means(start_pos_ind) = NaN;
    end


    %if any dieout cases, show these too
    if sum(which_sims_dieout)
        scatter(folder_inds(start_pos_ind)*ones(size(time_to_peak_all(which_sims_dieout))),...
            zeros(size(time_to_peak_all(which_sims_dieout))), ...
            50, dieout_colour, 'x', 'LineWidth', 2);
    end

end

%plot means
p1 = plot(folder_inds, peak_time_means, 'Marker', 'o', 'LineWidth', 1.5, 'Color', tpeak_colour);
p2 = plot(folder_inds, time_half_inv_means, 'Marker', 'o', 'LineWidth', 1.5, 'Color', t_05_colour);



%% formatting

ylim([0,plot_max])
xlim([0, gen_lims(end)])
xticks([0, gen_lims])
xlabel({'','Seed Cell Position'})
ylabel('Hours');


legend([p1, p2], {'$t_{\mathrm{peak}}$', '$t_{0.5}$'}, ...
    'Interpreter', 'latex', 'FontSize', round(0.9*fig_font_size), 'Location', 'SouthEast')


%optionally save
if save_stuff
    SaveAsPngAndFig(-1, 't_peak_and_t_05_seed_pos', fig_width, fig_ratio, fig_font_size);
    exportgraphics(gcf, 't_peak_and_t_05_seed_pos.png', 'Resolution', 300);
end


