function [ P_beta ] = prob_dist_power( P, beta, binom )
% Calculate the correctly normalized probability distribution when
% a population spike count distribution P(k) is raised to the power beta.

% BACKGROUND
% Given 3 neurons (say), the 2^3 possible states these can be in are:
% 000,001,010,100,011,101,110,111. If we take the neurons to be identical
% then we can say the probablity of seeing state 001 is the same as 010,
% and 100. It is then simpler to talk about population spike counts
% instead of particular states. In terms of population spike counts we 
% have 4 states: 0,1,2,3. If we denote by P(k) population spike count
% probability distribution, and P(x) the particular states prob dist then:
% P(k) = nchoosek*P(x) where |x|=k. i.e P(2) = (3choose2)*P(011).
% END BACKGROUND

% We first change from P(k) to P(x). We know that P(x) = P(k)/nchoosek. 
% We raise the probability distribution for this particular state to the 
% power beta.
P_pow_beta = (P./binom).^(beta);

% We must next correctly normalize this by dividing by the sum over all
% possible states. There are nchoosek possible ways that state x can occur
% so we account for all of these by multiplying by the correct binomial 
% coefficients in the sum.
Z_beta = sum(binom.*P_pow_beta);

% Finally we get P(k)^beta by multiplying P(x)^beta by the binomial 
% coefficients.
P_beta = binom.*P_pow_beta/Z_beta;

end

