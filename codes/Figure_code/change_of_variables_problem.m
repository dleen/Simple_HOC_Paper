%
% Probability change of variables problem
%

clear all
close all

% Values of lambda for linear filter (LF)
% and DG that produce the same statistics for each
lambda_LF = 0.22;
lambda_DG = 0.17;

% Time span over which to solve the equation
step = 0.01;
Tspan = -3:step:3;

% P_LF = probability distribution phi_LF(S) as given in
% the paper
% S_LF = the range of S, [-1, 1]
[P_LF, S_LF] = function_calc_P_LF_NL(lambda_LF, Tspan);

% P_DG = probability distribution phi_DG(s) in paper.
% pdf for the normal distribution
[P_DG, S_DG] = function_calc_P_DG(lambda_DG, Tspan);


% Initial condition
IC = 0;

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
LTildeS = 1-exp(-f_at_S_LF);


% The L(s) function for the DG. It's
% just the CDF:
gamma_dg = -1.28;
% rescale:
s_rescaled=(S_DG+gamma_dg)/sqrt(1-lambda_DG);
% Error function
DG_erf = erf(s_rescaled/sqrt(2));
% L(s) for DG:
DG = 0.5*(1+DG_erf);


%
% Plot results
%
%%
figure(3)
title('Comparing L(s) functions','fontsize',16)
% Plot support
p6=plot(S_DG,P_DG,'b--');
hold on
% Plot LF
p1=plot(S_LF,LTildeS,'r');
% Plot DG
plot(S_DG,DG,'b')
hold off
% %%
% % %%
% % % Plot comparison
% figure(1)
% % semilogy(S_DG, P_DG,'r')
% plot(S_DG, P_DG, 'r')
% hold on
% % Compare the change of variables
% % P_DG(s) = P_LF(f(s)) (df / ds)
% plot(S_LF, phi_LF_fs.*abs(dfds))
% hold off
% 
% figure(2)
% semilogy(S_DG, P_DG, 'r')
% hold on
% semilogy(S_LF, phi_LF_fs.*abs(dfds))
% % axis([-2 2 1e-3 1])
% hold off
% 
% figure(4)
% plot(rhs)

% L_function = [[c';ceif';clf';cising'],[N';N';N';N'],[ones(size(c'));2*ones(size(c'));3*ones(size(c'));4*ones(size(c'))]];
L_function = [[DG';LTildeS';P_DG'],[S_DG';S_LF';S_DG'],[ones(size(DG'));2*ones(size(DG'));3*ones(size(DG'))]];