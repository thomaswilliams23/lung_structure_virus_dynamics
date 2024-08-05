

%% tissue dimension parameters

%total depth of tree
total_tissue_size_x = 500;

%total circumference of tree
total_tissue_size_y = 64;

%total number of cells
all_cells_total = total_tissue_size_x*total_tissue_size_y;


%length of each tissue generation
generation_lengths = [100, 100, 100, 100, 100];

%number of tubes per generation
generation_num_tubes = 2.^(0:(length(generation_lengths)-1));

%circumferences of branch in each tissue generation
generation_widths =  round(total_tissue_size_y ./ generation_num_tubes);



%% model parameters

%gamma shape parameter for latent stage duration
latent_stages = 3;

%viral diffusion coefficient (CD^2\h^-1)
virus_diff = 100;

%rate of cell-to-cell infection (h^-1)
alpha_param = 1.839483;

%rate of extracellular viral infection ((TCID_50/ml)^-1 h^-1)
beta_param = 2.694176e-08;

%latently infected cell activation rate (h^-1)
gamma_param = 3.366934e-01;

%death rate of infected cells (h^-1)
delta_param = 8.256588e-02;

%extracellular virion production rate ((TCID_50/ml) h^-1)
p = exp(1.409457e+01);

%extracellular virion clearance rate (h^-1)
c = 4.313531e-01;


