r_range = 1:10000;

times = r_range*dt;

firing_rate = [[r_est(r_range);r_nl_full(r_range)],...
    [times';times'],...
    [ones(size(r_range))';2*ones(size(r_range))']];