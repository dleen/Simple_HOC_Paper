function [mu,rho,P_LF,cv,dbl_count] =...
    EIF_filter_dtp(gamma,lambda,sigma,N,alg,lnl,bin_size)

% Choice of algorithms:
% alg='maths for neuro'; % integrate firing rate until threshold.
% alg='mfn'; % independent implementation of above alg.
% alg='coin flip'; % r(t) \Delta t.
% alg='peche'; % Poisson.

%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
E0=gamma; % DC current
tau_m=5; % membrane time constant
tau_ref=3; % refractory period
deltat=3; % EIF "AP sharpness"
v_soft=-53; % EIF soft threshold
v_reset=-60; % Reset value
v_th=20; % EIF threshold
E1=sigma*sqrt(tau_m*lambda); % Perturbation current

sigmap=sigma/sqrt(2); % Value of sigma at which to...
% calculate the stationary firing rate
sigma_pert=sigma*sqrt(1-lambda)/sqrt(2); % Value of sigma at which to...
% calculate the linear filter.

Tmax = 200; % Max time lag for filter (ms)
dt = 0.1; % Integration time step (ms)
% Note: dtp < 0.02 gives NaN error

% Generate a vector of frequencies
dw = 1/2/Tmax; % Frequency step size
wmax = 1/2/dt; % Max frequency
w = -wmax:dw:(wmax-dw); % List of frequencies
%ind0 = find(abs(w) < 1e-6);
w(abs(w)<1e-6) = 1e-6;
bins = length(w);

%bin_size=10; % bin size in ms
num_of_bins=100000; % 

L=num_of_bins*bin_size/dt;

I_common = E1*randn(L,1)/sqrt(dt); % Common input

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

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

% Calculate the static non-linearity. See Ostoijic for details.
j=1;
E_range = (E0-10):0.051:(E0+10);
r0_func = zeros(size(E_range));
for i=E_range
    r0_func(j) = 1000*calc_Rate_cuEIF(i,sigmap,tau_ref,...
                    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    j=j+1;
end

r0_prime = gradient(r0_func,0.051);
r0_prime_E0 = interp1(E_range,r0_prime,E0,'spline')

n_max=400; 
D0 = sum(real(1000*A_t))*dt
nonlinear = zeros(n_max,1);
nonlinear_full = zeros(n_max,1);
E_span = linspace(-200,200,n_max); % 200Hz upper and lower...
                                   % limits for non-linearity.
for i=1:n_max
    nonlinear(i) = calc_Rate_cuEIF(E0+E_span(i)/D0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    nonlinear_full(i) = calc_Rate_cuEIF(E0+E_span(i)/r0_prime_E0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
end

% Calculate firing rates using the linear filter applied to common
% input. The non-linear firing rate is calculated by interpolating the
% static non-linearity.
r_est = r0 + conv(I_common,real(A_t),'same')*dt; % Linear firing rate.
r_nl = interp1(E_span/1000,nonlinear,r_est,'spline'); % LNL Cascade.
r_nl_full = interp1(E_span/1000,nonlinear_full,r_est,'spline'); % LNL Cascade.

% Generate spikes using the linear or non-linear firing rate:
if strcmp(lnl,'none')
    [mu,rho,P_LF,cv,dbl_count,~,~] =...
    Spikes_from_Lin_Filter(r_est,bin_size,N,alg,dt,num_of_bins);
elseif strcmp(lnl,'l')
    r_est=max(0,r_est);
    [mu,rho,P_LF,cv,dbl_count,~,~] =...
    Spikes_from_Lin_Filter(r_est,bin_size,N,alg,dt,num_of_bins);
elseif strcmp(lnl,'lnl')
    [mu,rho,P_LF,cv,dbl_count,~,~] =...
    Spikes_from_Lin_Filter(r_nl,bin_size,N,alg,dt,num_of_bins);
elseif strcmp(lnl,'lnl_full')
    [mu,rho,P_LF,cv,dbl_count,~,~] =...
    Spikes_from_Lin_Filter(r_nl_full,bin_size,N,alg,dt,num_of_bins);
end
end

