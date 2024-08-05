

%% time options

%time values (h)
dt = 0.01;
final_time = 24;

%early termination options
stop_when_all_infected = 0;
stop_if_dieout = 0;


%% plotting options

%plot cell grid during simulation?
vis_grid = 1;

%save out cell grid frames from simulation? (to make a movie of the
%simulation)
export_frames = 0;

%number of frames between visualisations of cell grid
frame_interval = round(1/dt);


%% number of cell lineages

num_lineages = 1;


%% boundary conditions at ends

use_toroidal_BCs = 0;


%% cell IDs

target_id = 1;
eclipse_ids = 2:(2+num_lineages-1);
infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_lineages - 1);
infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_lineages-1);
infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_lineages - 1);
infected_donor_dead_ids = (infected_donor_ids(end)+1):(infected_donor_ids(end)+1+num_lineages-1);
num_cell_types = infected_donor_dead_ids(end);


%set infected cell IDs as negative for faster lookup
infected_recip_ids = -infected_recip_ids;
infected_donor_ids = -infected_donor_ids;




