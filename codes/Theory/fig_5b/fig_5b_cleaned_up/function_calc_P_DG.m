function [P_DG, X] = function_calc_P_DG(lambda_DG)

X = -3:0.001:3;
P_DG = normpdf(X,0,sqrt(lambda_DG));