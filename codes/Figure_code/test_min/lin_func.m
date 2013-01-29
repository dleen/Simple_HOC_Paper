function [negative_log_likelihood] = lin_func( x, P_EIF, N, sigma )
%LIN_FUNC Summary of this function goes here
%   Detailed explanation goes here

[~,~,P_LF,~,~] = EIF_filter_test(x(1),x(2),sigma,N,'mfn','l',10,'normed');

negative_log_likelihood = -P_EIF'*log(P_LF);

end

