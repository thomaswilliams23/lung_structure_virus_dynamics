%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       LUNG STRUCTURE MODEL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements a CA virus model with two modes of viral spread

% Model structure:
%  - Cells: Target, Infected, Dead, arranged in hexagonal grid
%  - Extracellular virus: PDE

% This iteration is hard-coded to represent the spread of virus across a
% multiply-branched tube domain resembling the following: 
%
%
%                                        _____________
%                      ________________ /   ___G_2^1__)   ...
%                     /                :   /
%                    /      G_1^1      : ~ |              etc.
%    ---------------/     _____________:   \__________
%   (               :    /              \______G_2^2__)   ...
%   (    G_0^1      :~ ~|                _____________
%   (               :    \_____________ /   ___G_2^3__)   ...
%    ---------------\                  :   /
%                    \      G_1^2      : ~ |              etc.
%                     \________________:   \__________
%                                       \______G_2^4__)   ...
%
% Note that the branching is always on the right.


function [sim_out] = single_run_as_func(init_cells, out_times, params)


    addpath helper_functions


    %load default parameters
    biology_parameters;
    control_parameters;
    

    %read in and override default parameters:

    %control
    total_tissue_size_x = params.total_tissue_size_x;
    total_tissue_size_y = params.total_tissue_size_y;

    all_cells_total = total_tissue_size_x*total_tissue_size_y;

    generation_lengths = params.generation_lengths;
    
    final_time = params.final_time;
    stop_when_all_infected = params.stop_when_all_infected;
    stop_if_dieout = params.stop_if_dieout;
    
    vis_grid = params.vis_grid;
    frame_interval = params.frame_interval;
    export_frames = params.export_frames;

    num_lineages = params.num_lineages;

    use_toroidal_BCs = params.use_toroidal_BCs;
    
    
    
    %biological
    latent_stages = params.latent_stages;
    virus_diff = params.virus_diff;
    beta_param = params.beta_param;
    alpha_param = params.alpha_param;
    gamma_param = params.gamma_param;
    delta_param = params.delta_param;
    p = params.p;
    c = params.c;
    
    

    %set up output arrays
    prop_infected = 0*out_times;

    net_T = 0*out_times;
    net_E = 0*out_times;
    net_I = 0*out_times;
    net_V = 0*out_times;
    
    net_CC = 0*out_times;
    net_CF = 0*out_times;

    max_inv_depth = 0*out_times;
    mean_inv_depth = 0*out_times;

    if export_frames
        grid_frames = zeros(total_tissue_size_x, total_tissue_size_y, length(out_times));
        virus_frames = zeros(total_tissue_size_x, total_tissue_size_y, num_lineages, length(out_times));
        eclipse_wait_time_frames = zeros(total_tissue_size_x, total_tissue_size_y, length(out_times));
    end

    section_sums = zeros(length(generation_widths), length(out_times));
    
    %set up an output index
    curr_output_ind = 1;



    %set up cell IDs
    target_id = 1;
    eclipse_ids = 2:(2+num_lineages-1);
    infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_lineages - 1);
    infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_lineages-1);
    infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_lineages - 1);
    infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_lineages-1);
    num_cell_types = infected_donor_dead_ids(end);
    
    %infected IDs are negative for faster lookup
    infected_recip_ids = -infected_recip_ids;
    infected_donor_ids = -infected_donor_ids;


    %set up lineage order for visualisation
    lineage_order = randperm(num_lineages);

    %set up generation widths
    generation_num_tubes = 2.^(0:(length(generation_lengths)-1));
    generation_widths =  round(total_tissue_size_y ./ generation_num_tubes);
    generation_ends = [0, cumsum(generation_lengths)];

    %set up mesh for PDE
    mesh_parameters;
    make_diff_matrix;

    %initialise a count for infections from each mechanism
    tot_cf_inf = 0;
    tot_cc_inf = 0;



    %initialise, specify the initially infected cells
    cell_grid = target_id * ones(total_tissue_size_x, total_tissue_size_y);
    for cell_ind = 1:size(init_cells,1)
        cell_grid(init_cells(cell_ind,1),init_cells(cell_ind,2))=...
            infected_donor_ids(mod(cell_ind-1,num_lineages)+1);
    end

    %initialise eclipse stage durations
    eclipse_wait_times = 0*cell_grid;

    %initialise extracellular virus quantities
    if ~(virus_diff==Inf)
        v_ext=zeros(total_tissue_size_x, total_tissue_size_y, num_lineages);
    else
        v_ext=zeros(1,num_lineages);
    end



    %compute proportion of cells initially infected
    moi = sum(sum((cell_grid~=target_id)))/all_cells_total;



    %main loop
    t=0;
    time_last_infected = 0;
    stopping_criterion_reached=0;
    num_output_points = length(out_times);
    while t<=final_time && ~stopping_criterion_reached

        %finite diffusion
        if ~(virus_diff==Inf)
            %make next grid
            make_next_grid;
    
            %update v_ext
            update_v_ext;

        %infinite diffusion
        else
            %make next grid
            make_next_grid_global_virus;

            %update v_ext
            update_global_v_ext;
        end
        
        %update the grid
        cell_grid = next_grid;
        
        %update eclipse wait times
        eclipse_wait_times = eclipse_wait_times - dt;

        

        %visualise
        if vis_grid && ~mod(round(t/dt), frame_interval)
            pause(0.5)
            plot_grid;
        end


        %optional output at specified times
        if curr_output_ind <= length(out_times) &&...
                abs(t - out_times(curr_output_ind))<0.5*dt
    
            %save proportion of target cells
            net_T(curr_output_ind) = sum(sum((cell_grid == target_id)))/all_cells_total;
            
            %also virus
            if ~(virus_diff==Inf)
                net_V(curr_output_ind) = sum(sum(sum(v_ext)));
            else
                net_V(curr_output_ind) = sum(v_ext);
            end

            
            %loop over lineages
            for lineage = 1:num_lineages

                %save the proportion of recipient infected cells (incl. dead)
                prop_infected(curr_output_ind) = prop_infected(curr_output_ind) + ...
                    (sum(sum((cell_grid==infected_recip_ids(lineage)) + ...
                             (cell_grid==infected_recip_dead_ids(lineage)))))/((1-moi)*all_cells_total);
                    
                %save the proportion of target and infected cells
                net_E(curr_output_ind) = net_E(curr_output_ind) + sum(sum((cell_grid == eclipse_ids(lineage))))/all_cells_total;
                net_I(curr_output_ind) = net_I(curr_output_ind) + sum(sum((cell_grid == infected_recip_ids(lineage)) + ...
                    (cell_grid == infected_donor_ids(lineage))))/all_cells_total;


                %compute proportion infected per generation
                for gen = 1:length(generation_lengths)
                    section_sums(gen, curr_output_ind) = section_sums(gen, curr_output_ind) +...
                        sum(sum((cell_grid((generation_ends(gen)+1):generation_ends(gen+1),:) == infected_recip_ids(lineage)) + ...
                        (cell_grid((generation_ends(gen)+1):generation_ends(gen+1),:) == infected_donor_ids(lineage))))/all_cells_total;
                end
            end
            

            %optionally, save grid frames
            if export_frames
                grid_frames(:,:,curr_output_ind) = cell_grid;
                eclipse_wait_time_frames(:,:,curr_output_ind) = eclipse_wait_times;
                if ~(virus_diff==Inf)
                    virus_frames(:,:,curr_output_ind) = v_ext;
                end
            end
            
            %increment output index
            curr_output_ind = curr_output_ind+1;
        end
        
        

        %%%check if stopping criterion reached

        %case if all cells infected
        if stop_when_all_infected
            if sum(sum((cell_grid == target_id)))==0
                stopping_criterion_reached = 1;
                num_output_points=curr_output_ind-1;
                
                if (num_output_points < length(out_times))
                    disp('WARNING: did not reach all output times (all cells infected)')
                end
                
            end
        end

        %case if infection dies out
        if stop_if_dieout

            %update last time where there were any infected cells
            inf_cell_count = 0;
            for lineage = 1:num_lineages
                inf_cell_count = inf_cell_count + sum(sum((cell_grid == infected_recip_ids(lineage)) + ...
                    (cell_grid == infected_donor_ids(lineage))));
            end
            if inf_cell_count>0
                time_last_infected = t;
            end

            %check if we have spent enough time without any infected cells
            %to terminate
            TIME_FOR_DIEOUT = 10;
            if t-time_last_infected>TIME_FOR_DIEOUT
                stopping_criterion_reached = 1;
                num_output_points=curr_output_ind-1;
                
                if (num_output_points < length(out_times))
                    disp('WARNING: did not reach all output times (dieout)')
                end
            end
        end

            
        

        %increment time
        t=t+dt;
    end

    %compute proportion of infections from cell-to-cell infection
    if sum(sum((cell_grid == target_id)))==0
        pcc_this_sim = tot_cc_inf / (tot_cc_inf + tot_cf_inf);
    else
        pcc_this_sim = NaN;
    end
    

    %output information
    sim_out.prop_infected = prop_infected(1:num_output_points);
    sim_out.net_T = net_T(1:num_output_points);
    sim_out.net_E = net_E(1:num_output_points);
    sim_out.net_I = net_I(1:num_output_points);
    sim_out.net_V = net_V(1:num_output_points)/max(net_V);
    sim_out.section_sums = section_sums(:,1:num_output_points);
    sim_out.final_grid_state = abs(cell_grid);
    sim_out.num_output_points = num_output_points;
    sim_out.pcc_this_sim = pcc_this_sim;

    
    if export_frames
        sim_out.grid_frames = grid_frames;
        sim_out.eclipse_wait_time_frames = eclipse_wait_time_frames;
        if ~(virus_diff==Inf)
            sim_out.virus_frames = virus_frames;
        end
    end
    
end
    
    
    
    
    