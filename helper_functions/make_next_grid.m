


%updates the cell grid



next_grid = cell_grid;

%loop over cells
for i=1:total_tissue_size_x
    for j=1:total_tissue_size_y

        %check cell type:
        if cell_grid(i,j)==target_id  %cell is TARGET


            %%%compute probability of being infected on next time step:

            %initialise
            tot_neighbours = 6;
            tot_inf_neighbours = 0;

            %find locations of neighbours
            [i_neighbours, j_neighbours] = ...
                get_neighbour_indices(i,j, [total_tissue_size_x, total_tissue_size_y], ...
                generation_lengths, generation_widths);
            
            %calculate number of infected neighbours
            for neighbour = 1:tot_neighbours
                if ~isnan(i_neighbours(neighbour))
                    tot_inf_neighbours = tot_inf_neighbours + (cell_grid(i_neighbours(neighbour), j_neighbours(neighbour))<0);
                end
            end

            %and total virus at cell
            virus_at_cell = cell_area*v_ext(i,j,:);


            %overall prob of transition
            prob_t_to_e = 1 - exp(-(all_cells_total * beta_param * sum(virus_at_cell) + alpha_param*(tot_inf_neighbours/tot_neighbours))*dt);


            %compute probability of infection from either mechanism
            prob_CC_unscaled = 1-exp(-alpha_param*(tot_inf_neighbours/tot_neighbours)*dt);
            prob_CF_unscaled = 1-exp(-all_cells_total*beta_param*sum(virus_at_cell)*dt);
            if (prob_CC_unscaled || prob_CF_unscaled)
                prob_CC_inf = prob_t_to_e * (prob_CC_unscaled/(prob_CC_unscaled + prob_CF_unscaled));
                prob_CF_inf = prob_t_to_e * (prob_CF_unscaled/(prob_CC_unscaled + prob_CF_unscaled));
            else
                prob_CC_inf = 0;
                prob_CF_inf = 0;
            end

            %draw a random number
            rnd_draw = rand;

            
            %CELL-CELL INFECTION
            if (rnd_draw <= prob_CC_inf)

                %count a cc infection
                tot_cc_inf = tot_cc_inf + 1;

                %initialise
                inf_neighbours = zeros(1,num_lineages);

                %make a vector of num infected neighbours per lineage
                for neighbour = 1:tot_neighbours
                    if ~isnan(i_neighbours(neighbour))
                        lineage = get_infected_lineage(abs(cell_grid(i_neighbours(neighbour), j_neighbours(neighbour))), num_lineages);
                        if lineage>0
                            inf_neighbours(lineage) = inf_neighbours(lineage) + 1;
                        end
                    end
                end

                
                %quantify infection prob contributions from each lineage
                cumulative_unscaled_lineage_probs = zeros(1, num_lineages);
                cumulative_unscaled_prob = 0;
                for lineage = 1:num_lineages
                    cumulative_unscaled_prob = cumulative_unscaled_prob + ...
                        1 - exp(-(alpha_param*inf_neighbours(lineage)/tot_neighbours)*dt);
                    cumulative_unscaled_lineage_probs(lineage) = cumulative_unscaled_prob;
                end
                lineage_probs = prob_CC_inf * (cumulative_unscaled_lineage_probs/cumulative_unscaled_lineage_probs(end));


                %decide which lineage caused the infection
                for lineage = 1:num_lineages
                    if rnd_draw<=lineage_probs(lineage)
                        next_grid(i,j) = eclipse_ids(lineage);
                        wait_time = gamrnd(latent_stages, 1/(latent_stages*gamma_param));
                        eclipse_wait_times(i,j) = wait_time;
                        break
                    end
                end
                    
            
                
            %CELL-FREE INFECTION    
            elseif ((rnd_draw<=prob_t_to_e) && (rnd_draw>prob_CC_inf))

                %count a cf infection
                tot_cf_inf = tot_cf_inf + 1;

                %quantify infection prob contributions from each lineage
                cumulative_unscaled_lineage_probs = zeros(1, num_lineages);
                cumulative_unscaled_prob = 0;
                for lineage = 1:num_lineages
                    cumulative_unscaled_prob = cumulative_unscaled_prob + ...
                        1 - exp(-all_cells_total*beta_param*virus_at_cell(lineage)*dt);
                    cumulative_unscaled_lineage_probs(lineage) = cumulative_unscaled_prob;
                end
                lineage_probs = prob_CC_inf + prob_CF_inf * (cumulative_unscaled_lineage_probs/cumulative_unscaled_lineage_probs(end));    
                    
                
                %decide which lineage caused the infection
                for lineage = 1:num_lineages
                    if rnd_draw<=lineage_probs(lineage)
                        next_grid(i,j) = eclipse_ids(lineage);
                        wait_time = gamrnd(latent_stages, 1/(latent_stages*gamma_param));
                        eclipse_wait_times(i,j) = wait_time;
                        break
                    end
                end
                
                
            end
            

            
            
            
            
        elseif ismember(cell_grid(i,j),eclipse_ids)    %cell is ECLIPSE
            
            if eclipse_wait_times(i,j) < 0
                [~,lineage] = ismember(cell_grid(i,j), eclipse_ids);
                next_grid(i,j) = infected_recip_ids(lineage);
            end
            
            

        elseif ismember(cell_grid(i,j),infected_recip_ids) %cell is INFECTED (recipient)

            %compute prob of cell dying in next time step
            prob_i_to_d = 1-exp(-delta_param*dt);

            %draw random number, decide next state of cell
            if rand<prob_i_to_d %becomes dead
                [~,lineage] = ismember(cell_grid(i,j), infected_recip_ids);
                next_grid(i,j) = infected_recip_dead_ids(lineage);
            end
        
        
        
        elseif ismember(cell_grid(i,j),infected_donor_ids) %cell is INFECTED (donor)

            %compute prob of cell dying in next time step
            prob_i_to_d = 1-exp(-delta_param*dt);

            %draw random number, decide next state of cell
            if rand<prob_i_to_d %becomes dead
                [~,lineage] = ismember(cell_grid(i,j), infected_donor_ids);
                next_grid(i,j) = infected_donor_dead_ids(lineage);
            end
            
        end  %dead cells do nothing
        
    end
end