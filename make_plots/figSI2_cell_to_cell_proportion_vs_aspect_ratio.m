


%Makes SI Figure 3 from the Supplementary Information



%% data specs and initialisation

%first set text interpreter
set(0, 'DefaultTextInterpreter', 'latex')

%and add path
addpath plot_helpers




%data source - NOTE DIFFERENT TO THE ONE USED FOR FIGURE 2 BUT EQUIVALENT
source_folder = '../SIM_DATA_from_ms/SWEEP_OUTPUT_diffusion_and_aspect_ratio_sweep_with_Pcc/';
template_1 = '/tube_circ_';
template_2 = '_diff_';

%tube circumferences
circ_folders = {'4', '8', '16', '32', '64'};
circ_vals = 2*2.^(1:5);

%diffusion coefficients
diff_folders = {'0', '1', '10', '100', '1000', 'Inf'};
diff_legend = {'0', '1', '10', '100', '1000', '$\infty$'};

%tube lengths and aspect ratios
tot_cells = 4096;
tubes_lengths = tot_cells./circ_vals;
AR_vals = circ_vals./tubes_lengths;

%number iterations per AR-diffusion combination
num_reps = 10;

%colours
diff_plot_colours = min(1.3*summer(length(diff_folders)+3),1);

%save output?
save_stuff = 1;



%% main loop

%set up figure
plot_handles = zeros(length(diff_folders), 1);
figure


%loop over diffusion
for diff_ind = 1:length(diff_folders)

    %initialise
    Pcc_this_diff = zeros(length(circ_folders), 1);

    %loop over AR
    for circ_ind = 1:length(circ_folders)

        %load data
        load(strcat(source_folder, template_1, ...
            circ_folders{circ_ind}, template_2, ...
            diff_folders{diff_ind}, '/all_sim_data_this_fldr'));


        %compute mean Pcc (ignore NaNs)
        tot_Pcc = 0;
        num_non_nan = 0;
        for rep = 1:num_reps
            Pcc_this_rep = all_sim_data_this_fldr{rep};
            if ~isnan(Pcc_this_rep)
                tot_Pcc = tot_Pcc + Pcc_this_rep;
                num_non_nan = num_non_nan + 1;
            end
        end
        if num_non_nan>0
            Pcc_this_diff(circ_ind) = tot_Pcc/num_non_nan;
        else
            errmsg = sprintf("No valid Pcc counts for D=%s, circ=%s\n", ...
                diff_folders{diff_ind}, circ_folders{circ_ind});
            disp(errmsg)
            quit
        end


    end


    %add to the plot
    plot_handles(diff_ind) = semilogx(AR_vals, Pcc_this_diff,...
        'color', diff_plot_colours(diff_ind,:));
    set(plot_handles(diff_ind), 'linewidth', 1.5)
    set(plot_handles(diff_ind), 'marker', 'o')
    hold on

end


%formatting
lgd = legend(plot_handles, diff_legend, 'Interpreter', 'latex', ...
    'FontSize', 9, 'Location', 'EastOutside');
lgd.Title.String = '$D$';
lgd.Title.Interpreter = 'latex';
lgd.Title.FontSize = 9;
xlim([10^-2.6, 10^0.2])
ylim([0.84, 1.02])
xlabel('Aspect Ratio')
ylabel('Cell-to-cell infection proportion')
box off


%save
if save_stuff
    SaveAsPngAndFig([], 'Pcc_vs_AR', 10, 1.6, 9);
end



