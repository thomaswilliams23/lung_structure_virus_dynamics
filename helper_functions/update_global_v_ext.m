


%update the (global) extracellular virus



%compute a vector of infected proportion by viral lineage 
production_vector = 0*v_ext;

for lineage=1:num_lineages
    production_vector(lineage) = sum(sum( (cell_grid == infected_donor_ids(lineage)) + (cell_grid == infected_recip_ids(lineage))))/numel(cell_grid);
end
         

%update virus
v_ext = v_ext + dt*(p * production_vector - c*v_ext);