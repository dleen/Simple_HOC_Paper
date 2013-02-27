function [P_LF, X] = function_calc_P_LF_NL(lambda_LF, Tspan)
% 
%
%
%

%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
gamma=-60.31;

E0=gamma; % DC current
tau_m=5; % membrane time constant
tau_ref=3; % refractory period
deltat=3; % EIF "AP sharpness"
v_soft=-53; % EIF soft threshold
v_reset=-60; % Reset value
v_th=20; % EIF threshold
sigma=4.925; % variance
E1=sigma*sqrt(tau_m*lambda_LF); % Perturbation current

sigmap=sigma/sqrt(2); % Value of sigma at which to
% calculate the stationary firing rate
sigma_pert=sigma*sqrt(1-lambda_LF)/sqrt(2); % Value of sigma at which to
% calculate the linear filter.

Tmax = 200; % Max time lag for filter (ms)
dt = 0.1; % Integration time step (ms)
% Note: dt < 0.02 gives NaN error

% Generate a vector of frequencies
dw = 1/2/Tmax; % Frequency step size
wmax = 1/2/dt; % Max frequency
w = -wmax:dw:(wmax-dw); % List of frequencies
%ind0 = find(abs(w) < 1e-6);
w(abs(w)<1e-6) = 1e-6;
bins = length(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Richardson paper functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate steady state firing rate
r0=calc_Rate_cuEIF(E0,sigmap,tau_ref,...
    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);

% Calculate linear filter in frequency domain
Ahat=calc_Susc_cuEIF(w,E0,sigma_pert,tau_ref,...
    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01,r0);

% Fourier transform linear filter to temporal domain
[~,A_t]=inv_f_trans_on_vector(w,Ahat);
A_t=transpose(A_t);
A_t = fliplr(A_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Richardson functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the bins
% and the common input signal.
bin_size = 10; 
num_of_bins=200000; 
L=num_of_bins*bin_size/dt;
I_common = E1*randn(L,1)/sqrt(dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static nonlinearity. Ostoijic %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the static non-linearity. See Ostoijic for details.
j=1;
step_size=0.051;
E_range = (E0-10):step_size:(E0+10);
r0_func = zeros(size(E_range));
for i=E_range
    r0_func(j) = calc_Rate_cuEIF(i,sigmap,tau_ref,...
                    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    j=j+1;
end

r0_prime = gradient(r0_func,step_size);
r0_prime_E0 = interp1(E_range,r0_prime,E0,'spline');

n_max=400; 
D0 = sum(real(A_t))*dt;
nonlinear = zeros(n_max,1);
nonlinear_full = nonlinear;
E_span = linspace(-200/1000,200/1000,n_max); % 200Hz upper and lower...
                                   % limits for non-linearity.
for i=1:n_max
    % The difference between these two is the first uses
    % D0 in the expression: E0 + E_span/D0 where D0 is 
    % a constant calculated from the filter. See Ostojik
    % for more details.
    % The second one uses the more complicated form also found
    % in Ostojik. It uses r0_prime_E0 instead of D0 in the 
    % denominator. The difference is that r0_prime is calculated
    % using the steady state firing rate function above.
    nonlinear(i) = calc_Rate_cuEIF(E0+E_span(i)/D0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    nonlinear_full(i) = calc_Rate_cuEIF(E0+E_span(i)/r0_prime_E0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
end

%%%%%%%%%%%%%%%%%%%
% End of Ostoijic %
%%%%%%%%%%%%%%%%%%%

% Calculate firing rates using the linear filter applied to common
% input. The non-linear firing rate is calculated by interpolating the
% static non-linearity.
r_est = r0 + conv(I_common,real(A_t),'same')*dt; % Linear firing rate.
r_nl = interp1(E_span,nonlinear,r_est,'spline'); % LNL Cascade.
r_nl_full = interp1(E_span,nonlinear_full,r_est,'spline'); % LNL Cascade.

r_to_use = r_nl_full;


% Function from the paper
for i=1:(num_of_bins-floor(bin_size/dt))
    SS(i) =...
      sum(r_to_use(((i-1)*floor(bin_size/dt)+1):(i*floor(bin_size/dt))))*dt;
end


plf_s = -2:0.001:2;

% [P_LF, X] = ksdensity(SS,plf_s);
[P_LF, X] = ksdensity(SS, Tspan);
