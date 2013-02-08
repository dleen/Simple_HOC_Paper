%function [ dfds ] = prob_transf(s,f,lambda_DG,gm)
function [ dfds ] = prob_transf_back(s,f,plf_s,plf,lambda_DG )
% This function attempts to solve the change of variable ODE
% for a pdf i.e. http://en.wikipedia.org/wiki/
% Probability_density_function#Dependent_
% variables_and_change_of_variables
%
% i.e given df/ds = p1(s) / p2(f(s))
% solve for f

% Doesn't seem to be working quite right though.

    %prob_dg = interp1(pdg_s,pdg,s,'spline');
    prob_lf = interp1(plf_s,plf,f,'spline');
    
    prob_dg = normpdf(-0.5-s,0,sqrt(lambda_DG));
    
    % Return derivative
    dfds = prob_dg/prob_lf;
end