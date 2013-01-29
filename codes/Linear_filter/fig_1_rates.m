r_est_none = r0 + conv(I_common,real(A_t),'same')*dt;

plot_length = 1:1000;

time_axis = cumsum(plot_length*dt);

r_none_plot = r_est_none(plot_length);
r_plot = r_est(plot_length);


firing_rate_R = [[r_plot;r_none_plot],[time_axis';time_axis'],...
    [ones(size(plot_length'));2*ones(size(plot_length'))]];