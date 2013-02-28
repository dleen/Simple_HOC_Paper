function [ mu, rho, sigma, P_DG ] = DG_statistics( lambda, gamma, N )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate DG Spikes.
% Columns of 'DG_spike_stack' represent the N neurons and rows represent
% the 'time' steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_samples = 500000;

% Following Macke notation:
DG_spike_stack = sign(sqrt(1-lambda)*randn(num_samples,N) + ...
                 sqrt(lambda)*randn(num_samples,1)*ones(1,N) + gamma);
             
% We'll use 0's and 1's as it makes the statistics etc much easier.
DG_spike_stack(DG_spike_stack == -1) = 0; 

% Temp vector to hold the number of occurences.
occurr = zeros(N+1,1);
% Sum across the rows to get the population spike count.
holder = sum(DG_spike_stack,2); 

% In this loop we count the number of times each spike count occurs.
% E.g. If holder = (1,1,1,4,3,3,7,2)^t.
% ind = find(holder == 3) will give the indices where the number 3 can 
% be found: ind = (5,6) and by taking the length of ind we get the 
% number of times 3 occurs. length(ind) = 2. This is the number of
% occurrences of the spike count 3.
for i=0:N
    tempind = find(holder == i); % ^^^^
    occurr(i+1) = length(tempind); % count the number of occurrences
end

% Finally, return the population spike count probability. This P(k) gives
% the probability that a state with k spikes will occur.
P_DG = occurr/sum(occurr);

% Statistics of the DG.
mu = mean(mean(DG_spike_stack)); % mean firing rate.
sigma = cov(DG_spike_stack); % find the covariance of the output.
%rho = mean(sigma(1,2:N))/mean(diag(sigma)); % defined in the paper.
rho=((sum(sum(sigma))-sum(diag(sigma)))/(N*(N-1)))/mean(diag(sigma));

end

