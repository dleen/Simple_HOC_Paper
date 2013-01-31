% 
%
%
%

clear all;
close all;

lambda_DG = 0.17;
%lambda_DG = 0.242;

%lambda_LF = 0.215;
%lambda_LF = 0.37;
lambda_LF = 0.22;

gamma_dg = -1.28;


%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
gamma=-60.31;
%gamma=-60.78;
%gamma=-59.05;
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

%%
warning('off')

bin_size = 10; 

%%%%%%
%%%%%%
%%%%%%
num_of_bins=200000; 
L=num_of_bins*bin_size/dt;
I_common = E1*randn(L,1)/sqrt(dt);

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

% Calculate firing rates using the linear filter applied to common
% input. The non-linear firing rate is calculated by interpolating the
% static non-linearity.
r_est = r0 + conv(I_common,real(A_t),'same')*dt; % Linear firing rate.
r_nl = interp1(E_span,nonlinear,r_est,'spline'); % LNL Cascade.
r_nl_full = interp1(E_span,nonlinear_full,r_est,'spline'); % LNL Cascade.

length(r_est(r_est<0))/length(r_est)

%%
close all

r_temp = r_nl_full;

for i=1:(num_of_bins-floor(bin_size/dt))
    SS(i) =...
      sum(r_temp(((i-1)*floor(bin_size/dt)+1):(i*floor(bin_size/dt))))*dt;
end

s_dg = -3:0.001:3;
pdg = normpdf(s_dg,0,sqrt(lambda_DG));

s_lf = s_dg;

plf_s = -1:0.001:1;

[plf,plf_s] = ksdensity(SS,plf_s);


%plf_s = plf_s';
%plf = abs(plf);

ind = plf_s<0;
plf(ind) = 0;
%ind = plf_s>0.5;
%plf(plf<1e-4) = 1e-4;

%gm = gamfit(SS);
%plf_gam=gampdf(plf_s(ind),gm(1),gm(2));
%plf(ind) = plf_gam;

pdg_s = s_dg; 

xx=(s_dg+gamma_dg)/sqrt(1-lambda_DG);
erftest = erf(xx/sqrt(2));
DG = 0.5*(1+erftest);

%Tspan = [-50:0.01:3]; % Solve from t=1 to t=5
%Tspan = [-3:0.01:3];
Tspan = -0.5:0.001:1;

ind = find(abs(s_dg+0.5)<1e-4);

IC = DG(ind); % y(t=0) = 1

%[t ff] = ode45(@(s1,f)prob_transf(s1,f,pdg_s,pdg,plf_s,plf),Tspan,IC); % Solve ODE
options = odeset('NonNegative',[],'RelTol',1e-3);
%[t ff] = ode23(@(s1,f)prob_transf(s1,f,lambda_DG,gm),Tspan,IC); % Solve ODE
[t ff] = ode23(@(s1,f)prob_transf(s1,f,plf_s,plf,lambda_DG),Tspan,IC); % Solve ODE

Tspan = -1:0.001:-0.5;
%[t1 ff1] = ode23(@(s1,f)prob_transf(s1,f,lambda_DG,gm),Tspan,IC); % Solve ODE
[t1 ff1] = ode23(@(s1,f)prob_transf_back(s1,f,plf_s,plf,lambda_DG),Tspan,IC); % Solve ODE

t  = [t1(1:(end-1));t];
ff = [ff1(1:(end-1));ff];

ffs = interp1(t,ff,plf_s,'spline');
phifs = interp1(plf_s,plf,ffs,'spline');
%dfds = interp1(t(1:(end)),gradient(ff(1:end),t(1:end)),plf_s1,'spline');
dfds = interp1(t,gradient(ff,t),plf_s,'spline');

% % DG
% s_dg1 = -3:0.01:3;
% pdg1 = normpdf(s_dg,0,sqrt(lambda_DG));

% 
EIF = 1-exp(-plf_s);
EIF_transf = 1-exp(-ffs);
% EIF_plot = EIF;
EIF = max(0,EIF);
% 
% xx=(ffs+gamma_dg)/sqrt(1-lambda_DG);
% erftest = erf(xx/sqrt(2));
% DG_nl = 0.5*(1+erftest);
% 
figure(2)
title('Comparing L(s) function','fontsize',16)
p6=plot(s_dg,DG);
hold on
p1=plot(plf_s,EIF_transf,'r');
plot(pdg_s,pdg,'b--')

figure(11)
plot(s_dg,pdg,'r')
hold on
plot(plf_s,phifs.*dfds)

figure(12)
title('Comparing L(s) function','fontsize',16)
p6=semilogy(s_dg,DG);
hold on
p1=semilogy(plf_s,EIF_transf,'r');
semilogy(pdg_s,pdg,'b--')
axis([-3 3 1e-7 1])

figure(14)
plot(t,ff)
title('Solution of ODE','fontsize',16)

figure(15)
plot(t,log(ff))
title('Solution of ODE','fontsize',16)

figure(6)
plot(plf_s,plf)
title('Probability density of nonlin S','fontsize',16)

N=100;

P_DG_func = @(kk) nchoosek(N,kk)*...
     sum(pdg.*((1-DG).^(N-kk)).*((DG).^kk))*(s_dg(2)-s_dg(1));
P_LF_func = @(kk) nchoosek(N,kk)*...
     sum(plf.*((1-EIF).^(N-kk)).*((EIF).^kk))*(plf_s(2)-plf_s(1));


 for i=0:N
     P_DG_list(i+1) = P_DG_func(i);
     P_LF_list(i+1) = P_LF_func(i);
 end

figure(10)
plot(0:N,P_DG_list,'b')
hold on
plot(0:N,P_LF_list,'r')

figure(3)
pp1=semilogy(0:N,P_DG_list,'b'); 
hold on
 pp2=semilogy(0:N,P_LF_list,'r');
axis([0 N 1e-4 1])
xlabel('Population Spike Count  k','fontsize',16)
ylabel('P(k)','fontsize',16)
title(['Comparison of Linear Filter simulation with theory, with DG',...
    10,'Mean firing rate \mu = 0.1, correlation \rho = 0'],...
    'fontsize',16)

figure(1)
semilogy(s_dg,pdg,'r')
hold on
semilogy(plf_s,phifs.*dfds)
axis([-2 2 1e-3 1])