function [negative_log_likelihood,gradient_of_negative_log_likelihood] =...
    neg_log_like_binom(...
    param_list,feature_matrix,mean_feature_d,prob_d,binom)

prob_m_unnormalized = binom.*exp(feature_matrix * param_list);  
%this is the prob. of each state occurring under max-ent distrib.
%Col vector with length = number of states
prob_m = prob_m_unnormalized/sum(prob_m_unnormalized) ;

mean_feature_m =  prob_m' * feature_matrix ; %row vector

%negative_log_likelihood = - sum (prob_d .* log( prob_m)) ;
negative_log_likelihood = -prob_d'*log(prob_m);

gradient_of_negative_log_likelihood =...
    -( mean_feature_d - mean_feature_m )';
%not obvious -- comes from form of max ent distribs






