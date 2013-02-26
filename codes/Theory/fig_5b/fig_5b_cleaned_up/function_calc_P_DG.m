function [P_DG, X] = function_calc_P_DG(lambda_DG, Tspan)

% X = -3:0.001:3;
X = Tspan;
P_DG = normpdf(X,0,sqrt(lambda_DG));