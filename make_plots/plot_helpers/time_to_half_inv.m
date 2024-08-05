function time_val = time_to_half_inv(data_vec)

    %find m such that half the area under data_vec is under data_vec(1:m)

    data_cum_sum = cumsum(data_vec)/sum(data_vec);
    time_val = find(data_cum_sum>0.5, 1, 'first');

end