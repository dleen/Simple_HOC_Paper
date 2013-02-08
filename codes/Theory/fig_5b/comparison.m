clear all; close all; clf;

%lambda_DG = 0.132;
%lambda_DG = 0;
lambda_DG = 0.242;
%lambda_LF = 0.248;
lambda_LF = 0.532;
%lambda_LF = 0.532;
%gamma_dg = 0;
gamma_dg = -1.28;

%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
gamma=-61.1;
%gamma=-61.1;
E0=gamma; % DC current
%tau_m=5; % membrane time constant
tau_m=5; % membrane time constant
tau_ref=3; % refractory period
deltat=3; % EIF "AP sharpness"
v_soft=-53; % EIF soft threshold
v_reset=-60; % Reset value
v_th=20; % EIF threshold
sigma=6.295;
%sigma=4.5;
E1=sigma*sqrt(tau_m*lambda_LF); % Perturbation current

sigmap=sigma/sqrt(2); % Value of sigma at which to...
% calculate the stationary firing rate
sigma_pert=sigma*sqrt(1-lambda_LF)/sqrt(2); % Value of sigma at which to...
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




%% Non linearity

j=1;
step_size = 0.006; % Weird NaN bug at E = -56 mV.
E_range = (E0-10):step_size:(E0+10);
r0_func = zeros(size(E_range));
for i=E_range
    r0_func(j) = calc_Rate_cuEIF(i,sigmap,tau_ref,...
                    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    j=j+1;
end

r0_prime = gradient(r0_func,step_size);
r0_prime_E0 = interp1(E_range,r0_prime,E0,'spline');


D0 = sum(real(A_t))*dt;

n_max=400;
nonlinear = zeros(n_max,1);
nonlinear_full = nonlinear;
E_span = linspace(-200/1000,200/1000,n_max);
for i=1:n_max
    nonlinear_full(i) = calc_Rate_cuEIF(E0+E_span(i)/r0_prime_E0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    nonlinear(i) = calc_Rate_cuEIF(E0+E_span(i)/D0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
end

%%
clf


s = -4:0.005:8;

xx=(s+gamma_dg)/sqrt(1-lambda_DG);
erftest = erf(xx/sqrt(2));


DG = 0.5*(1+erftest);

deltacapt = 10; % 100 ms bin size is the best?

bbeta = sum(real(A_t))*dt;

%beta = sum(real(A_t))*dt*deltacapt;
%gamma_EIF = r0*deltacapt;

% alpha = sigma^2*tau_m*lambda_LF;
% CC = sqrt(alpha/lambda_DG)*bbeta;

alpha = sigma^2*tau_m;
CC = sqrt(alpha)*bbeta %CC = 0.003
r0
ss = -4:0.005:8;

% r_nl =...
%     interp1(E_span,nonlinear,r0+bbeta*sqrt(alpha/lambda_DG)*ss,'spline');
% 
% r_nl_full =...
%     interp1(E_span,nonlinear_full,r0+bbeta*sqrt(alpha/lambda_DG)*ss,'spline');
r_est = r0+CC*ss;
r_nl =...
    interp1(E_span,nonlinear,r0+CC*ss,'spline');

r_nl_full =...
    interp1(E_span,nonlinear_full,r0+CC*ss,'spline');


%EIF = 1-exp(-beta*ss-gamma_EIF);
%EIF_lin = gamma_EIF + beta*ss;

EIF_lin  = 1-exp(-r_est*deltacapt);
EIF      = 1-exp(-r_nl*deltacapt);
EIF_full = 1-exp(-r_nl_full*deltacapt);

figure(1)
plot(s,DG);
hold on
plot(ss,EIF_lin,'k')
plot(ss,EIF,'g')
plot(ss,EIF_full,'r')
axis([-4 8 -0.1 1.1])
legend('DG','Linear','Nonlinear 1','Nonlinear 2','location','northwest')

LDG0 = interp1(s,DG,0);
L0 = interp1(s,EIF_lin,0);
LNL0 = interp1(s,EIF,0);

N=100;
P_DG = @(k) nchoosek(N,k)*((1-LDG0)^(N-k))*(LDG0)^k;

P_LF_lin = @(k) nchoosek(N,k)*((1-L0)^(N-k))*(L0)^k;

for i=0:100
    P_DG_list(i+1) = P_DG(i);
end

figure(2)
plot(0:100,P_DG_list)
