Tspan = [-0.1 1]; % Solve from t=1 to t=5
IC = 0; % y(t=0) = 1
[t ff] = ode45(@(s,f)prob_transf(s,f,pdg_s,pdg,plf_s,plf),Tspan,IC); % Solve ODE