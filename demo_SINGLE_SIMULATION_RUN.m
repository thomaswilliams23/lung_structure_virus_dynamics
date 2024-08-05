%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       LUNG STRUCTURE MODEL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements a CA virus model with two modes of viral spread

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



addpath helper_functions

%load default params
prms = make_default_param_struct;



%%% CHANGE PARAMETERS HERE

%visualising
prms.vis_grid = 0;
prms.frame_interval = 10;
prms.export_frames = 0;

%number of lineages
prms.num_lineages = 1;

%timing
prms.final_time = 300;
prms.stop_when_all_infected = 0;
prms.stop_if_dieout = 1;

%dimensions
prms.total_tissue_size_x = 500;
prms.total_tissue_size_y = 64;
prms.all_cells_total = prms.total_tissue_size_x * prms.total_tissue_size_y;
prms.generation_lengths = [100, 100, 100, 100, 100];
assert(sum(prms.generation_lengths) == prms.total_tissue_size_x);

%biology
prms.latent_stages = 3;
prms.virus_diff = 100;
prms.alpha_param = 1.83948353e+00;
prms.beta_param = 2.69417636e-08;
prms.gamma_param = 3.366934e-01; 
prms.delta_param = 8.256588e-02;
prms.p = exp(1.409457e+01);
prms.c = 4.313531e-01;
%%%




%specify any data output times
vis_dt = 0.5;
out_times = 0:vis_dt:(prms.final_time-vis_dt);


%specify initial cells
init_cells = [1, randi(prms.total_tissue_size_y)];

%run simulation
[sim_out] = single_run_as_func(init_cells, out_times, prms);



