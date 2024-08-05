


%given a cell's id, retrieves its lineage number



function species = get_infected_lineage(id, num_lineages)

    %set up ids
    eclipse_ids = 2:(2+num_lineages-1);
    infected_recip_ids = (eclipse_ids(end) + 1):(eclipse_ids(end) + 1 + num_lineages - 1);
    infected_recip_dead_ids = (infected_recip_ids(end)+1):(infected_recip_ids(end)+1 + num_lineages-1);
    infected_donor_ids = (infected_recip_dead_ids(end) + 1):(infected_recip_dead_ids(end) + 1 + num_lineages - 1);

    %find if the id belongs to one of the infected ids
    [~,s1] = ismember(id, infected_recip_ids);
    [~,s2] = ismember(id, infected_donor_ids);

    species = s1+s2;
end
