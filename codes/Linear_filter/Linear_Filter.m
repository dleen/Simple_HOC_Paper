%clear all; clf; close all;
tic

%%%%%%%%%%%%%%%%%%%%%%%%
% Compile if necessary %
%%%%%%%%%%%%%%%%%%%%%%%%
%mex -largeArrayDims calc_Rate_cuEIF.cpp
%mex -largeArrayDims calc_Susc_cuEIF.cpp

%%%%%%%%%%%%%%%%%%%%%%%%
% Linear filter values %
%%%%%%%%%%%%%%%%%%%%%%%%
Tmax = 200; % Maximum time lag (ms)
dt = 0.1; % Integration time step (ms)
% Note: dt < 0.02 gives NaN error

% Generate a vector of frequencies
dw = 1/2/Tmax;
wmax = 1/2/dt;
w = -wmax:dw:(wmax-dw);
ind0 = find(abs(w) < 1e-6);
w(ind0) = 1e-6;
bins = length(w);

%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
sigma=4.925;
gamma_dc=-60.31;
E0=gamma_dc;
tau_m=5;
tau_ref=3;
deltat=3;
v_soft=-53;
v_reset=-60;
v_th=20;

lambda=0.215;

N=100; % Number of neurons

bin_size=10; % bin size in ms
num_of_bins=100000; % 
L=num_of_bins*bin_size/dt;

% Common input signal
E1=sigma*sqrt(tau_m*lambda);
I_common = E1*randn(L,1)/sqrt(dt);

%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithms available %
%%%%%%%%%%%%%%%%%%%%%%%%
%algorithm='maths for neuro'; 
algorithm='mfn';
%algorithm='coin flip';
%algorithm='peche';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Linear Filter in frequency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We can either use lambda dependent term in front of the independent
% input noise, realising that this will lead to the mean firing rate
% being dependent on lambda:
sigma_pert=sigma*sqrt((1-lambda)/2);
%sigma_pert=sigma*sqrt(1-lambda);
% Or we can set the term independent of lambda (the sqrt(2) is necessary
% to cancel a surplus sqrt(2) in Richardson's equation:
sigmap=sigma/sqrt(2);
%sigmap=sigma;

% Calculate stationary firing rate
r0=calc_Rate_cuEIF(E0,sigmap,tau_ref,...
    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);

% Calculate fourier transformed filter
Ahat=calc_Susc_cuEIF(w,E0,sigma_pert,tau_ref,...
    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01,r0);

% Fourier transform back
[time_ms,A_t]=inv_f_trans_on_vector(w,Ahat);
A_t=transpose(A_t);
A_t = fliplr(A_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linearity calculation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second Ostojic method
j=1;
step_size = 0.011;
E_range = (E0-5):step_size:(E0+5);
r0_func = zeros(size(E_range));
for i=E_range
    r0_func(j) = calc_Rate_cuEIF(i,sigmap,tau_ref,...
                    v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    j=j+1;
end
r0_prime = gradient(r0_func,step_size);
% Constant needed for second Ostojic method
r0_prime_E0 = interp1(E_range,r0_prime,E0,'spline');

% First Ostojic method constant
D0 = sum(real(A_t))*dt;

n_max=400;
nonlinear = zeros(n_max,1);
nonlinear_full = nonlinear;
E_span = linspace(-100/1000,200/1000,n_max); %range of nonlinearity in kHz
% Calculate non-linearity
for i=1:n_max
    % calculate nonlinearity using second Ostojic method
    nonlinear_full(i) = calc_Rate_cuEIF(E0+E_span(i)/r0_prime_E0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
    % calculate nonlinearity using first Ostojic method
    nonlinear(i) = calc_Rate_cuEIF(E0+E_span(i)/D0,...
        sigmap,tau_ref,v_reset,v_soft,v_th,deltat,tau_m,-100,0.01);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate spikes and calc statistics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear estimate of firing rate
r_est = r0 + conv(I_common,real(A_t),'same')*dt;
%r_est = max(0,r_est);

% Linear-Nonlinear estimate of firing rate using first Ostojic
r_nl = interp1(E_span,nonlinear,r_est,'spline');

% Linear-Nonlinear estimate of firing rate using second Ostojic
r_nl_full = interp1(E_span,nonlinear_full,r_est,'spline');

% Calculate spikes using whichever firing rate
% Change first input to r_est or r_nl etc
[mu,rho,P_Th,cv,dbl,t_trunc,n] =...
    Spikes_from_Lin_Filter(r_nl_full,bin_size,N,...
    algorithm,dt,num_of_bins);

%%
%%%%%%%%%%%%%%%%%
% Print results %
%%%%%%%%%%%%%%%%%
disp(algorithm)
disp('Steady state firing rate (Hz)')
disp(1000*r0)
disp('Mean firing rate')
disp(mu)
disp('Corr coeff')
disp(rho)
disp('CV')
disp(cv)
disp('Double count')
disp(dbl)
disp('Prob dist')
disp(P_Th)

%%
%%%%%%%%%
% Plots %
%%%%%%%%%
figure(1)
subplot(3,2,1)
plot(w((bins/2-bins/4):(bins/2+bins/4)),...
    abs(Ahat((bins/2-bins/4):(bins/2+bins/4))))
xlabel('Frequency (kHz)')
ylabel('$|\hat{R}(\omega)|$ (kHz/mV)','Interpreter','latex')
axis('tight')

subplot(3,2,2)
plot(time_ms((bins/2-bins/4):(bins/2+bins/4)),...
    1000*real(A_t((bins/2-bins/4):(bins/2+bins/4))))
xlabel('Time (ms)')
ylabel('R(t) Linear Filter (ms)')
axis('tight')

subplot(3,2,[5 6])
plot((0:100)*dt,1000*r_nl_full(1:101),'b')
xlabel('Time (ms)')
hold on
plot((0:100)*dt,1000*r_nl(1:101),'g')
plot((0:100)*dt,1000*r_est(1:101),'r')
plot(0:3*dt:100*dt,0,'r.','MarkerSize',0.5)
hold off
ylabel('Firing rate (Hz)')
axis('auto')

subplot(3,2,[3 4])
plot((0:3000)*dt,1000*r_nl_full(1:3001),'b')
hold on
plot((0:1000)*dt,1000*r_nl(1:1001),'g')
plot((0:3000)*dt,1000*r_est(1:3001),'r')
xlabel('Time (ms)')
hold on
plot(0:3*dt:1000*dt,0,'r.','MarkerSize',0.5)
hold off
ylabel('Firing rate (Hz)')
axis('auto')

figure(2)
plot(0:N,P_Th)
title('N = 10 neurons')
xlabel('Population spike count')
ylabel('P(n)')
legend('Linear Filter')

figure(3)
plot(1000*E_span,1000*nonlinear_full)
hold on
plot(1000*E_span,1000*nonlinear,'g')
plot(-100:200,(-100:200)+...
1000*interp1(E_span,nonlinear,0,'spline'),'r--')
legend('Full Non-linearity','Non-linearity',...
    'Unit slope','location','northwest')
axis('auto')
title('Non-linearity')
xlabel('Input Frequency')
ylabel('Output Frequency')
hold off

%%%%%%%
% End %
%%%%%%%
toc