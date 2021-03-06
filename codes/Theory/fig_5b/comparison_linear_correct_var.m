clear all; close all;

lambda_DG = 0.132;
%lambda_DG = 0.242;
%lambda_DG = 0.19;
%lambda_DG = 0.17;
%lambda_DG = 0.16;
%lambda_DG = 0.157;
lambda_LF = 0.25;
%lambda_LF = 0.575;
%lambda_LF = 0.532;
%lambda_LF = 0.365;
%lambda_LF = 0.37;
%lambda_LF = 0.45;
%lambda_LF = 0.64;
gamma_dg = -1.28;
%gamma_dg = -0.84;
%gamma_dg = -0.251;
%gamma_dg = 0;

%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
gamma=-60.3;

%gamma=-52.43;
%gamma=-54.5;
%gamma=-56.36;
%gamma=-58.3;
%gamma=-61.1;
%gamma=-61.45;
E0=gamma; % DC current
tau_m=5; % membrane time constant
tau_ref=3; % refractory period
deltat=3; % EIF "AP sharpness"
v_soft=-53; % EIF soft threshold
v_reset=-60; % Reset value
v_th=20; % EIF threshold
%sigma=6.295;
sigma=6.19;
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

%%
%clf
warning('off')

s = -10:0.005:10;

% DG
xx=(s+gamma_dg)/sqrt(1-lambda_DG);
erftest = erf(xx/sqrt(2));
DG = 0.5*(1+erftest);

%bin_size = 10; 
bin_size = 10; 

k=1;
for i=1:(bins-floor(bin_size/dt))
    B(k) = sum(real(...
    A_t((i):(i+floor(bin_size/dt)))... %fix the range here
        ))*dt;
    k=k+1;
end

var_2 = sum(B.^2)*dt;
E1^2*var_2
%%%%%%
%%%%%%
%%%%%%
num_of_bins=100000; 
L=num_of_bins*bin_size/dt;
I_common = E1*randn(L,1)/sqrt(dt);
r_est = conv(I_common,real(A_t),'same')*dt;

for i=1:(num_of_bins-floor(bin_size/dt))
    SS(i) =...
        sum(r_est(((i-1)*floor(bin_size/dt)+1):(i*floor(bin_size/dt))))*dt;
end
var(SS)
k=1;
for i=0:(bins/2-1)
    B_test(k) = sum(real(...
    A_t((bins/2-i):(bins/2+floor(bin_size/dt)-i))... %fix the range here
        ))*dt;
    k=k+1;
end
var_2_test = sum(B_test.^2)*dt;
%%%%%%
%%%%%%
%%%%%%

if (lambda_DG == 0 || lambda_LF == 0)
    alpha = sigma^2*tau_m*var_2;
    CC = sqrt(alpha); %CC = 0.003
else
    alpha = E1^2*var_2;
    CC = sqrt(alpha/lambda_DG);
end

EIF = 1-exp(-r0*bin_size-CC*s);
EIF_plot = EIF;
EIF = max(0,EIF);
%%
figure(1)
p6=plot(s,DG);
hold on
p1=plot(s,EIF,'r');
p2=plot(s,EIF_plot,'r--');

axis([-5 8 -0.1 1.1]);
p4=plot(0,-0.1:0.01:1.1);
p5=plot(-5:0.05:8,0.1,'k');
p3=plot(s,normpdf(s,0,sqrt(lambda_DG)),'--');
legend([p6 p1 p3],'DG','Linear est',...
    'PDF    N(0,\lambda)','location','northwest')
title(['Probability of spike conditioned on common input',...
    10,'Mean firing rate \mu = 0.1, correlation \rho = 0'],...
    'fontsize',16)
xlabel('s','fontsize',16)
ylabel('L(s) = P(spike|s)','fontsize',16)

% Values at s = 0;
LDG0 = interp1(s,DG,0,'spline');
L0 = interp1(s,EIF,0,'spline');

pdftest = normpdf(s,0,sqrt(lambda_DG));

N=64;

if (lambda_DG == 0 || lambda_LF == 0)
    P_DG_func = @(kk) nchoosek(N,kk)*...
        ((1-LDG0).^(N-kk)).*((LDG0).^kk);
    P_LF_func = @(kk) nchoosek(N,kk)*...
        ((1-L0).^(N-kk)).*((L0).^kk);
else
    P_DG_func = @(kk) nchoosek(N,kk)*...
        sum(pdftest.*((1-DG).^(N-kk)).*((DG).^kk))*0.005;
    P_LF_func = @(kk) nchoosek(N,kk)*...
        sum(pdftest.*((1-EIF).^(N-kk)).*((EIF).^kk))*0.005;
end

for i=0:N
    P_DG_list(i+1) = P_DG_func(i);
    P_LF_list(i+1) = P_LF_func(i);
end

figure(2)
plot(0:N,P_DG_list,'b')
hold on
plot(0:N,P_LF_list,'r')

figure(3)
pp1=semilogy(0:N,P_DG_list,'b');
hold on
pp2=semilogy(0:N,P_LF_list,'r');
axis([0 N 1e-5 1])
xlabel('Population Spike Count  k','fontsize',16)
ylabel('P(k)','fontsize',16)
title(['Comparison of Linear Filter simulation with theory, with DG',...
    10,'Mean firing rate \mu = 0.1, correlation \rho = 0'],...
    'fontsize',16)