%
% Probability change of variables problem
%

clear all
close all

% Values of lambda for linear filter (LF)
% and DG that produce the same statistics for each
lambda_LF = 0.22;
lambda_DG = 0.17;

% P_LF = probability distribution phi_LF(S) as given in
% the paper
% S_LF = the range of S, [-1, 1]
[P_LF, S_LF] = function_calc_P_LF_NL(lambda_LF);
% P_DG = probability distribution phi_DG(s) in paper.
% pdf for the normal distribution
[P_DG, S_DG] = function_calc_P_DG(lambda_DG);

% Time span over which to solve the equation
Tspan = -0.5:0.001:1;
% Initial condition
IC = 0;

% Solve the initial condition using ode45
% We are solving the ODE:
% df/ds = phi_DG(s) / phi_LF(f(s))
[s, f] = ode45(@(s,fs)right_hand_side_ratio(s, fs, P_LF, S_LF,... 
    P_DG, S_DG), Tspan, IC);


% Interpolate the values of f at the points S_LF
f_at_S_LF = interp1(s, f, S_LF,'spline');
% Interpolate the values of P_LF at f(s)
phi_LF_fs = interp1(S_LF, P_LF, f_at_S_LF,'spline');
% Interpolate the gradient of f at S_LF
dfds = interp1(s, gradient(f, s), S_LF,'spline');

% Plot comparison
figure(1)
semilogy(S_DG, P_DG,'r')
hold on
% Compare the change of variables
% P_DG(s) = P_LF(f(s)) (df / ds)
semilogy(S_LF, phi_LF_fs.*abs(dfds))
axis([-2 2 1e-3 1])