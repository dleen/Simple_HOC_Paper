%function [ dfds ] = prob_transf(s,f,lambda_DG,gm)
function [ dfds ] = prob_transf(s,f,plf_s,plf,lambda_DG )
    %prob_dg = interp1(pdg_s,pdg,s,'spline');
    prob_lf = interp1(plf_s,plf,f,'spline');
    
    %prob_lf = gampdf(f,gm(1),gm(2));
    prob_dg = normpdf(s,0,sqrt(lambda_DG));
    
    dfds = prob_dg/prob_lf;
end