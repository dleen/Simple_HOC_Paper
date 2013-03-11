%
% Probability change of variables problem
%

clear all
%close all

% Values of lambda for linear filter (LF)
% and DG that produce the same statistics for each
lambda_LF = 0.22;
lambda_DG = 0.17;

% Time span over which to solve the equation
step = 0.01;
Tspan = -3:step:3;

% P_LF = probability distribution phi_LF(S) as given in
% the paper
% S_LF = the range of S, [-1, 1]
[P_LF, S_LF] = function_calc_P_LF_NL(lambda_LF, Tspan);

% P_DG = probability distribution phi_DG(s) in paper.
% pdf for the normal distribution
[P_DG, S_DG] = function_calc_P_DG(lambda_DG, Tspan);


% i = 1
% for ic = -0.01:0.001:0
%     best_ic_vals(i) = find_best_IC_func(ic, P_LF, S_LF, P_DG, S_DG, Tspan, step);
%     i = i + 1;
% end

hold on
ic = 0;
i = 1;
curr = -1;
while ~isnan(curr)
    curr = find_best_IC_func(ic, P_LF, S_LF, P_DG, S_DG, Tspan, step);
    best_ic_vals(i) = curr;
    ic_vals(i) = ic;
    i = i + 1;
    ic = ic - 0.001;
end

best_ic_vals(isnan(best_ic_vals)) = 0;
best_ic_vals(isinf(best_ic_vals)) = 0;
best_ic_vals(best_ic_vals > 1000) = 0;

plot(ic_vals(1:15), best_ic_vals(1:15), 'b')
