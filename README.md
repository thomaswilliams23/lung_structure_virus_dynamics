Code used for Williams et al., "Accounting for the geometry of the respiratory tract in viral infections"

Simulation data which can be used with this code can be found on Figshare at DOI:10.26188/26499811


All code is written in MATLAB. This repository is structured as follows. 

In the main folder:
- "demo_SINGLE_SIMULATION_RUN" runs a demonstration simulation of the model. Parameters can be changed within the script where indicated to alter the geometry and dynamics of the simulation.
- Scripts with the prefix "RUN_SIMULATIONS_" run parameter sweeps on the model (exact scope of the scripts are specified within the script preambles) and can be used to generate the kind of datasets used in the model.
- "single_run_as_func" calls the model as a function (necessary for the other scripts in this folder).

In "make_plots" we include code to generate all of the figures in the manuscript (figure numbers indicated in the file name). These scripts call data which can be generated using the scripts mentioned above; alternatively, the datasets we used for the manuscript are available to download from FigShare (doi:10.26188/26499811). Scripts in subfolder "plot_helpers" are used for figure generation, of these, "patchline", "SaveAsPngAndFig" and "violin" are third-party scripts freely available from MATLAB File Exchange, credit is given in these scripts.

In "helper_functions" we include additional scripts called by the main function, "single_run_as_func". Of these, "circles" is a third-party script freely available from MATLAB File Exchange, credit is given in this script.
