


% given cell indices i,j, returns a list of i and j indices of the six
% neighbours of cell i, j (NaN if on an edge)
% order is [BL, TL, B, T, BR, TR]



function [i_inds, j_inds] = get_neighbour_indices(i,j,sheet_dims, ...
    generation_lengths, generation_circs)


%unpack dimensions
sheet_size_x = sheet_dims(1);
sheet_size_y = sheet_dims(2);

%number of neighbours is fixed
num_neighbours = 6;

%initialise
i_inds = zeros(num_neighbours,1);
j_inds = zeros(num_neighbours,1);


%% deal with default first, cell in body of tissue (will fix branch cases later)
if mod(i,2) %odd column, bottom-aligned

    if i==1
        i_inds([1,2]) = NaN;
        j_inds([1,2]) = NaN;
    else
        %BL
        i_inds(1) = i-1;
        j_inds(1) = mod(j-2,sheet_size_y)+1;

        %TL
        i_inds(2) = i-1;
        j_inds(2) = j;
    end

    %B
    i_inds(3) = i;
    j_inds(3) = mod(j-2,sheet_size_y)+1;

    %T
    i_inds(4) = i;
    j_inds(4) = mod(j,sheet_size_y)+1;

    if i==sheet_size_x
        i_inds([5,6]) = NaN;
        j_inds([5,6]) = NaN;
    else
        %BR
        i_inds(5) = i+1;
        j_inds(5) = mod(j-2,sheet_size_y)+1;

        %TR
        i_inds(6) = i+1;
        j_inds(6) = j;
    end


else %even column, top-aligned

    
    if i==1
        i_inds([1,2]) = NaN;
        j_inds([1,2]) = NaN;
    else
        %BL
        i_inds(1) = i-1;
        j_inds(1) = j;

        %TL
        i_inds(2) = i-1;
        j_inds(2) = mod(j,sheet_size_y)+1;
    end

    %B
    i_inds(3) = i;
    j_inds(3) = mod(j-2,sheet_size_y)+1;

    %T
    i_inds(4) = i;
    j_inds(4) = mod(j,sheet_size_y)+1;

    if i==sheet_size_x
        i_inds([5,6]) = NaN;
        j_inds([5,6]) = NaN;
    else
        %BR
        i_inds(5) = i+1;
        j_inds(5) = j;

        %TR
        i_inds(6) = i+1;
        j_inds(6) = mod(j,sheet_size_y)+1;
    end
end


%% now fix edges


%simple case: single tube
if length(generation_lengths)==1
    return

%harder case: branching tubes
else

    %calculate which generation the cell is in
    gen=1;
    while i > sum(generation_lengths(1:gen))
        gen=gen+1;
    end
    
    
    %easy case, on generation with no split
    if gen==1
        return
    
    %harder case, on generation with splits
    else
    
        %tube circumference for this generation
        tube_circ = generation_circs(gen);
    
        %check if on a split
        if mod(j,tube_circ)>1
            return
    
        else
            %at top of tube
            if mod(j,tube_circ)==0
    
                if mod(i,2) %odd column, bottom-aligned
                    j_inds(4) = mod(j_inds(4)-tube_circ-1,sheet_size_y) + 1; %T
                else %even column, top-aligned
                    j_inds(2) = mod(j_inds(1)-tube_circ-1,sheet_size_y) + 1; %TL
                    j_inds(4) = mod(j_inds(4)-tube_circ-1,sheet_size_y) + 1; %T
                    j_inds(6) = mod(j_inds(5)-tube_circ-1,sheet_size_y) + 1; %TR
                end
    
            %at bottom of tube
            elseif mod(j,tube_circ)==1
    
                if mod(i,2) %odd column, bottom-aligned
                    j_inds(1) = mod(j_inds(1)+tube_circ-1,sheet_size_y) + 1; %BL
                    j_inds(3) = mod(j_inds(3)+tube_circ-1,sheet_size_y) + 1; %B
                    j_inds(5) = mod(j_inds(5)+tube_circ-1,sheet_size_y) + 1; %BR
                else %even column, bottom-aligned
                    j_inds(3) = mod(j_inds(3)+tube_circ-1,sheet_size_y) + 1; %B
                end
    
            end
        end
    end
    
    
    %now just tidy up any NaNs
    for k = 1:num_neighbours
        if isnan(i_inds(k))
            j_inds(k) = NaN;
        end
    end
    
    
end

end

