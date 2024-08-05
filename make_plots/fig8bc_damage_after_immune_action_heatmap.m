


%Makes Figure 8(b), (c) and (d) from the manuscript and SI Figure 3 from
%the Supplementary Information



%% setup

%first set text intepreter
set(0, 'DefaultTextInterpreter', 'latex')

%and path
addpath plot_helpers


%data sources
branching_data_open_edge = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_top_and_bottom_seed/open_edge';
branching_data_branching_edge = '../SIM_DATA_from_ms/SWEEP_OUTPUT_branching_top_and_bottom_seed/branching_edge';
tube_data_wide = '../SIM_DATA_from_ms/SWEEP_OUTPUT_wide_and_narrow_tube_fixed_length/tube_circ_64';
tube_data_narrow = '../SIM_DATA_from_ms/SWEEP_OUTPUT_wide_and_narrow_tube_fixed_length/tube_circ_4';


%data specs
num_reps = 100;


%contour levels
highlight_contour = 0.5;
contours_for_comparison = 0:0.1:1;


%sweep ranges
thresh_cutoffs = 0.01*(1:50);
act_time_sweep = 1*(0:80);

%legends
contour_comparison_legend = {'Source Open Edge', 'Source Branching Edge'};
Finf_05_comparison_legend = {'Source Open Edge', 'Source Branching Edge', 'Gen. 1 Tube', 'Gen. 5 Tube'};

%time
max_t = 300;
dt = 0.1;


%colours
open_edge_colour = [217, 95, 2]/256; 
branching_edge_colour = [102, 166, 30]/256;


%save?
save_stuff = 1;



%% analyse data


%loop over branching case and non branching case
for branching_version = [0,1]
    for starts_on_left_edge = [0,1]

        %initialise
        all_final_damage = zeros(length(thresh_cutoffs), length(act_time_sweep), num_reps);

        %load in data for tube simulations (these are packaged in one file)
        if ~branching_version
            if starts_on_left_edge
                load(strcat(tube_data_wide, '/all_sim_data_this_fldr'));
            else
                load(strcat(tube_data_narrow, '/all_sim_data_this_fldr'));
            end
        end


        %loop over reps
        for rep = 1:num_reps

            %load data
            if branching_version
                if starts_on_left_edge
                    load(strcat(branching_data_open_edge, '/sim_data_', num2str(rep)));
                    this_sim = sim_data_this_rep;
                else
                    load(strcat(branching_data_branching_edge, '/sim_data_', num2str(rep)));
                    this_sim = sim_data_this_rep;
                end
            else
                this_sim = all_sim_data_this_fldr{rep};
            end



            %loop over thresholds
            for thresh_ind = 1:length(thresh_cutoffs)

                %set threshold
                thresh_cutoff = thresh_cutoffs(thresh_ind);

                %loop over activation times
                for act_time_ind = 1:length(act_time_sweep)

                    %set activation time
                    act_time = act_time_sweep(act_time_ind);


                    %find if and when the cumulative infected proportion passes
                    %this threshold
                    [thresh_exceeded,ind_of_thresh] = max((this_sim.prop_infected>thresh_cutoff));


                    %decide if a dieout case
                    DIEOUT_CUTOFF = 0.9;
                    not_dieout = (this_sim.prop_infected(end)>DIEOUT_CUTOFF);


                    %if a valid case, compute total damage after activation time
                    if thresh_exceeded && not_dieout

                        %find time of threshold and index of immune activation
                        time_to_thresh = ind_of_thresh*dt;
                        immune_activation_ind = round(time_to_thresh/dt) + ...
                            round(act_time/dt);

                        %case where immune activation happens in time
                        if immune_activation_ind <= this_sim.num_output_points
                            all_final_damage(thresh_ind, act_time_ind, rep) = ...
                                1-this_sim.net_T(immune_activation_ind);

                        %case where doesn't
                        else
                            all_final_damage(thresh_ind, act_time_ind, rep) = 1;
                        end


                    %if invalid case (dieout or threshold not reached), write NaN
                    else

                        all_final_damage(thresh_ind, act_time_ind, rep) = NaN;
                    end

                end
            end
        end


        %now rename the array
        if branching_version
            if starts_on_left_edge
                all_final_damage_open = all_final_damage;
            else
                all_final_damage_branching = all_final_damage;
            end
        else
            if starts_on_left_edge
                all_final_damage_wide_tube = all_final_damage;
            else
                all_final_damage_narrow_tube = all_final_damage;
            end
        end

    end
end


%% plot heatmaps


%%%open edge:

%take a mean
final_damage_means_open = squeeze(mean(all_final_damage_open, 3, 'omitnan'));

%plot contours
figure
h = contourf(act_time_sweep, thresh_cutoffs, final_damage_means_open, 'EdgeColor','none');
hold on

%plot the highlight contour (0.5)
[Ch, Hh] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_open, ...
    highlight_contour*[1,1], 'Color', open_edge_colour, 'LineWidth', 1.5);
clabel(Ch, Hh, 'manual', 'FontSize', 9, 'Color', open_edge_colour, 'Interpreter', 'latex')

%heatmap colours
colours = colormap('sky');
colormap(0.3 + 0.7*colours)
colorbar

%formatting
xlabel('$t_{\mathrm{act}}$ (h)');
ylabel('$F^{\mathrm{thresh}}$');
title('Source Open Edge', 'interpreter', 'latex')


%optionally save
if save_stuff
    SaveAsPngAndFig([], 'final_damage_heatmap_open_edge', 7, 4/3, 9);
end



%%%branching edge:

%take a mean
final_damage_means_branching = squeeze(mean(all_final_damage_branching, 3, 'omitnan'));

%plot contours
figure
h = contourf(act_time_sweep, thresh_cutoffs, final_damage_means_branching, 'EdgeColor','none');
hold on

%plot the highlight contour (0.5)
[Ch, Hh] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_branching, ...
    highlight_contour*[1,1], 'Color', branching_edge_colour, 'LineWidth', 1.5);
clabel(Ch, Hh, 'manual', 'FontSize', 9, 'Color', branching_edge_colour, 'Interpreter', 'latex')

%heatmap colours
colours = colormap('sky');
colormap(0.3 + 0.7*colours)
colorbar

%formatting
xlabel('$t_{\mathrm{act}}$ (h)');
ylabel('$F^{\mathrm{thresh}}$');
title('Source Branched Edge', 'interpreter', 'latex')


%optionally save
if save_stuff
    SaveAsPngAndFig([], 'final_damage_heatmap_branched_edge', 7, 4/3, 9);
end




%% compare all contours for the open and branched cases

%plot
figure
[c1, h1] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_open, contours_for_comparison, ...
    'Color', open_edge_colour, 'LineWidth', 1.5);
clabel(c1, h1, 'Color', open_edge_colour);
hold on
[c2, h2] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_branching, contours_for_comparison, ...
    'Color', branching_edge_colour, 'LineWidth', 1.5);
clabel(c2, h2, 'Color', branching_edge_colour);


%dummy lines for legend
p1 = plot([0,0], [0,0], 'Color', open_edge_colour);
p2 = plot([0,0], [0,0], 'Color', branching_edge_colour);

%formatting
legend([p1, p2], contour_comparison_legend, 'FontSize', 9, 'Interpreter', 'latex')
title('$F^{\infty}$ Contours')
xlabel('$t_{\mathrm{act}}$ (h)')
ylabel('$F^{\mathrm{thresh}}$')

%optionally save
if save_stuff
    SaveAsPngAndFig([], 'all_contours_comparison', 10, 1.15, 10);
end



%% now compare the highlight contour (0.5)

%plot
figure
hold on
[cm1, ch1] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_open, ...
    highlight_contour * [1,1], 'Color', open_edge_colour, 'LineWidth', 2);
[cm2, ch2] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_branching, ...
    highlight_contour * [1,1], 'Color', branching_edge_colour, 'LineWidth', 2);


%%% shade regions

%extract coords of highlight contours
contour_open_edge_x = cm1(1,2:end);
contour_open_edge_y = cm1(2,2:end);

contour_branching_edge_x = cm2(1,2:end);
contour_branching_edge_y = cm2(2,2:end);

%convert to splines
open_edge_spline = spline(contour_open_edge_x, [contour_open_edge_y(1:end-1), 0]);
branching_edge_spline = spline(contour_branching_edge_x, contour_branching_edge_y);

%find the x values for the contour intersection and the axis intercept
cand_x_vals = linspace(0, act_time_sweep(end), 200);
[~,crossover_ind] = min(abs(ppval(open_edge_spline,cand_x_vals)-ppval(branching_edge_spline, cand_x_vals)));
[~,x_intercept_ind_1] = min(abs(ppval(open_edge_spline,cand_x_vals)));

%region 1: y<open_contour && y<branching_contour
bottom_region_x = cand_x_vals(1:x_intercept_ind_1);
patch('XData', [bottom_region_x, fliplr(bottom_region_x)], 'YData', ...
    [min(ppval(open_edge_spline, bottom_region_x), ppval(branching_edge_spline, bottom_region_x)),...
    0*bottom_region_x], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4)


%region 2: y<open_contour && y>branching_contour
left_region_x = cand_x_vals(1:crossover_ind);
patch('XData', [left_region_x, fliplr(left_region_x)], 'YData', ...
    [ppval(open_edge_spline, left_region_x), fliplr(ppval(branching_edge_spline, left_region_x))],...
    'FaceColor', open_edge_colour, 'EdgeColor', 'none', 'FaceAlpha', 0.2)


%region 3: y>open_contour && y<branching_contour
right_region_x = cand_x_vals(crossover_ind:end);
patch('XData', [right_region_x, fliplr(right_region_x)], 'YData', ...
    [max(ppval(open_edge_spline, right_region_x),0*right_region_x), ...
    fliplr(ppval(branching_edge_spline, right_region_x))],...
    'FaceColor', branching_edge_colour, 'EdgeColor', 'none', 'FaceAlpha', 0.2)


%region 4: y>open_contour && y>branching_contour
top_region_x = cand_x_vals;
patch('XData', [top_region_x, fliplr(top_region_x)], 'YData', ...
    [max(ppval(open_edge_spline, top_region_x), ppval(branching_edge_spline, top_region_x)),...
    0*top_region_x+0.5], 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.4)


%annotate the regions
text(12, 0.35, 'I', 'Interpreter', 'latex', 'FontSize', 12)
text(60, 0.05, 'II', 'Interpreter', 'latex', 'FontSize', 12)



%now process the tube data
final_damage_means_wide_tube = squeeze(mean(all_final_damage_wide_tube, 3, 'omitnan'));
final_damage_means_narrow_tube = squeeze(mean(all_final_damage_wide_tube, 3, 'omitnan'));

%and add to the plot
[cm3, ch3] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_wide_tube, ...
    highlight_contour * [1,1], 'Color', open_edge_colour, 'LineWidth', 1.5, 'LineStyle', '--');
[cm4, ch4] = contour(act_time_sweep, thresh_cutoffs, final_damage_means_narrow_tube, ...
    highlight_contour * [1,1], 'Color', branching_edge_colour, 'LineWidth', 1.5, 'LineStyle', '--');

%formatting
legend([ch1, ch2, ch3, ch4], Finf_05_comparison_legend, 'Interpreter','latex')
xlabel('$t_{\mathrm{act}}$ (h)');
ylabel('$F^{\mathrm{thresh}}$');
title('$F^{\infty}=0.5$ Contour Comparison')
box on



%optionally save
if save_stuff
    SaveAsPngAndFig([], 'Finf_05_contour_comparison', 10, 1.15, 10)
    exportgraphics(gcf, 'Finf_05_contour_comparison.png', 'Resolution',300)
end
