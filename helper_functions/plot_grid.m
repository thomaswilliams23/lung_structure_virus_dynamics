


%code for plotting the state of the grid



%% colour scheme

%lineage colours
lineage_colours = [127,201,127; 190,174,212; 253,192,134; 255,255,153]/256;

%target and dead cell colours, depending on whether the simulation is
%finished or ongoing
if exist('plotting_final_grid', 'var')
    if plotting_final_grid
        target_colour = [0.8,0.8,0.8];
    else
        target_colour = [0.8, 0.8, 0.8];
        dead_colour = [0.3,0.3,0.3];
    end
else
    plotting_final_grid = 0;
    target_colour = [0.8, 0.8, 0.8];
    dead_colour = [0.3,0.3,0.3];
end


%% preprocess the grid to mark with lineage number
if plotting_final_grid
    get_spec_num = @(id) get_infected_lineage_incl_dead(id, num_lineages);
else
    get_spec_num = @(id) get_infected_lineage(id, num_lineages);
end
cell_grid_lineages = arrayfun(get_spec_num, abs(cell_grid));


%% make the figure

%initialise - depends if we want to destroy the figure or just update
if plotting_final_grid
    figure
else
    figure(1)
end

%loop over lineage (assumes everything is infected)
for lineage = 1:num_lineages

    %find cells with this cell type
    [x_vals, y_vals] = find(cell_grid_lineages==lineage);

    %offset for hexagonal structure
    y_vals = y_vals + 0.5*(mod(x_vals,2));

    %adjust x_vals for true hexagonal packing
    x_vals = (sqrt(3)/2) * x_vals;

    %if any cells with this cell type
    if ~isempty(x_vals)

        %plot circular cells in the correct loaction
        circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', lineage_colours(lineage,:), ...
                'EdgeColor', 'none');
        hold on
    end

end

%% now plot target and dead cells
if plotting_final_grid
    %set up target id
    target_id = 0;

    %plot target cells
    [x_vals, y_vals] = find(cell_grid_lineages==target_id);
    y_vals = y_vals + 0.5*(mod(x_vals,2));
    x_vals = (sqrt(3)/2) * x_vals;
    if ~isempty(x_vals)
        circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', target_colour, ...
                'EdgeColor', 'none');
    end

else
    %set up cell ids
    target_id = 1;
    eclipse_ids = 2:(2+num_lineages-1);
    infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_lineages - 1);
    infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_lineages-1);
    infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_lineages - 1);
    infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_lineages-1);
    num_cell_types = infected_donor_dead_ids(end);

    %plot target cells
    [x_vals, y_vals] = find(abs(cell_grid)<infected_recip_ids(1));
    y_vals = y_vals + 0.5*(mod(x_vals,2));
    x_vals = (sqrt(3)/2) * x_vals;
    if ~isempty(x_vals)
        circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', target_colour, ...
                'EdgeColor', 'none');
    end

    %plot dead cells
    which_cells_dead = 0*cell_grid;
    for lineage = 1:num_lineages
        which_cells_dead = which_cells_dead + (abs(cell_grid)==infected_donor_dead_ids(lineage)) + ...
                                              (abs(cell_grid)==infected_recip_dead_ids(lineage));
    end
    [x_vals, y_vals] = find(which_cells_dead);
    y_vals = y_vals + 0.5*(mod(x_vals,2));
    x_vals = (sqrt(3)/2) * x_vals;
    if ~isempty(x_vals)
        circles(x_vals, y_vals, 0.5*ones(size(x_vals)), 'FaceColor', dead_colour, ...
                'EdgeColor', 'none');
    end
end



%% plot branch divisions
generation_starts = cumsum(generation_lengths);
num_branches = 2.^(0:(length(generation_lengths)-1));
for gen = 1:(length(generation_starts)-1)
    x_vals_on_edge_unscaled = (generation_starts(gen):generation_starts(gen+1));
    x_vals_on_edge = sqrt(3)/2 * (x_vals_on_edge_unscaled+0.5);

    for split_num = 1:(num_branches(gen+1)-1) %this number of splits
        y_vals_on_edge = split_num*size(cell_grid,2)/num_branches(gen+1) ...
            + 0.5 + 0.5*(mod(x_vals_on_edge_unscaled,2));
        plot(x_vals_on_edge, y_vals_on_edge, 'k', 'LineWidth', 1.2)
    end
end


%% formatting


xticks([])
yticks([])

xlim([0,size(cell_grid,1)+0.5])
ylim([0,size(cell_grid,2)+1])

axis off

x_size = ((sqrt(3)/2)*size(cell_grid,1));
y_size = size(cell_grid,2);
ax = gca;
ax.PlotBoxAspectRatio = [1,y_size/x_size,1];

fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [0.1, 0.1, 0.8, 0.8];

hold off

