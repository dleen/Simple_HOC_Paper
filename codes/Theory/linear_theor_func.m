function [P_LF] =...
    linear_theor_func(gamma,lambda,sigma,N,bin_size)

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
% Note: dt < 0.02 gives NaN error

% Generate a vector of frequencies
dw = 1/2/Tmax; % Frequency step size
wmax = 1/2/dt; % Max frequency
w = -wmax:dw:(wmax-dw); % List of frequencies
%ind0 = find(abs(w) < 1e-6);
w(abs(w)<1e-6) = 1e-6;
bins = length(w);

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

warning('off')

s = -10:0.005:10;

k=1;
for i=1:(bins-floor(bin_size/dt))
    B(k) = sum(real(...
    A_t((i):(i+floor(bin_size/dt)))... %fix the range here
        ))*dt;
    k=k+1;
end

var_2 = sum(B.^2)*dt;

if (lambda == 0)
    alpha = sigma^2*tau_m*var_2;
    CC = sqrt(alpha); %CC = 0.003
else
    alpha = E1^2*var_2;
    CC = sqrt(alpha);
    %CC = sqrt(alpha/lambda_DG);
end

%EIF = 1-exp(-r0*bin_size-CC*s);
EIF = 1-exp(-r0*bin_size-s);
EIF_plot = EIF;
EIF = max(0,EIF);


% Values at s = 0;
L0 = interp1(s,EIF,0,'spline');

pdftest = normpdf(s,0,sqrt(alpha));

if (lambda == 0)
    P_LF_func = @(kk) nchoosek(N,kk)*...
        ((1-L0).^(N-kk)).*((L0).^kk);
else
    P_LF_func = @(kk) nchoosek(N,kk)*...
        sum(pdftest.*((1-EIF).^(N-kk)).*((EIF).^kk))*0.005;
end

for i=0:100
    P_LF(i+1) = P_LF_func(i);
end

end
