function [ nm ] = find_best_IC_func( IC, P_LF, S_LF, P_DG, S_DG, Tspan, step )


% Initial condition
% IC = 0;

% Solve the initial condition using ode45
% We are solving the ODE:
% df/ds = phi_DG(s) / phi_LF(f(s))
% [s, f] = ode45(@(s,f)right_hand_side_ratio(s, f, P_LF, S_LF,... 
%     P_DG, S_DG), Tspan, IC);

[s, f, rhs] = my_euler(P_LF, S_LF, P_DG, S_DG, Tspan, step, IC);


% Interpolate the values of f at the points S_LF
f_at_S_LF = interp1(s, f, S_LF,'spline');
% Interpolate the values of P_LF at f(s)
phi_LF_fs = interp1(S_LF, P_LF, f_at_S_LF,'spline');
% Interpolate the gradient of f at S_LF
dfds = interp1(s, gradient(f, s), S_LF,'spline');


%
% Calculate L(s) functions
%

% The L^\tilde(S) function for the Linear Filter:
% LTildeS = 1-exp(-f_at_S_LF);


% % The L(s) function for the DG. It's
% % just the CDF:
% gamma_dg = -1.28;
% % rescale:
% s_rescaled=(S_DG+gamma_dg)/sqrt(1-lambda_DG);
% % Error function
% DG_erf = erf(s_rescaled/sqrt(2));
% % L(s) for DG:
% DG = 0.5*(1+DG_erf);

% norm(P_DG)
% norm(f)
% norm(phi_LF_fs.*abs(dfds))

nm = norm(phi_LF_fs.*abs(dfds) - P_DG);

end