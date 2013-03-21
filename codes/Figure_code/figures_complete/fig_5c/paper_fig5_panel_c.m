clear all
close all

lambda_LF = 0.37;

%%%%%%%%%%%%%%
% LIF values %
%%%%%%%%%%%%%%
gamma=-60.78;
E0=gamma; % DC current
tau_m=5; % membrane time constant
tau_ref=3; % refractory period
deltat=3; % EIF "AP sharpness"
v_soft=-53; % EIF soft threshold
v_reset=-60; % Reset value
v_th=20; % EIF threshold
sigma=4.925;
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

warning('off')

bin_size = 10; 

%%%%%%
num_of_bins=100000; 
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

r_temp = r_nl_full;

for i=1:(num_of_bins-floor(bin_size/dt))
    SS(i) =...
    sum(r_temp(((i-1)*floor(bin_size/dt)+1):(i*floor(bin_size/dt))))*dt;
end

s_dg = -3:0.01:5;
s_lf = s_dg;

[plf,plf_s] = ksdensity(SS,s_lf);

EIF = 1-exp(-plf_s);
EIF = max(0,EIF);

P_LF_func = @(kk,N) nchoosek(N,kk)*...
    sum(plf.*((1-EIF).^(N-kk)).*((EIF).^kk))*(plf_s(2)-plf_s(1));

% Number of neurons.
N = 4:8:100;

% Temperature
T = 1;

for j=1:length(N)
    % Calculate the correct binomial coefficients.
    binom = zeros(N(j)+1,1);
    for i=0:N(j)
        binom(i+1) = nchoosek(N(j),i);
    end

    Q = zeros(N(j),1);
    for i=0:N(j)
        Q(i+1) = P_LF_func(i,N(j));
    end
    
    % Calculate the DG probability distribution.
    [mu,rho,~,P_DG] = DG_statistics(0.242,-1.28,N(j));
    P = P_DG;
    
    M{1} = importdata(strcat('../numerical_data/fig2mu/fig_2mu_',...
        int2str(N(j)),'_',num2str(0.1),'.dat'),' ');
    % PP contains the probability distributions.
    % Each row represents a single run at mean mu and corr rho.
    PP = M{1}(:,3:N(j)+3);
    mu_1 = M{1}(:,N(j)+4);
    rho_1 = M{1}(:,N(j)+5);

    P_EIF = mean(PP,1)';
    
    states = [(0:N(j))',(0:N(j))'.^2];
    mean_feature = P_EIF'*states;
            param_list_init = [0;0];

        % MATLAB Optimization Toolbox minimization function.
        % Uses Optimization toolbox unconstrained minimization function.
        options = optimset('GradObj','on','LargeScale','on',...
            'Display','off',...
            'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-8,'TolX',1e-8);
    %     options = optimset('GradObj','on',...
    %         'Display','final-detailed',...
    %     'MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-10,'TolX',1e-10);
        % We set GradObj = on as we supply the gradient.
        % We set LargeScale = on as that is the algorithm that uses the
        % user supplied gradient.
        % Display just shows some accuracy output.
        % The final options are tolerances etc.

        % The minimization function
        [param_list,fval,exitflag,output] =...
              fminunc(@(x)neg_log_like_binom...
              (x,states,mean_feature,Q,binom),param_list_init,options);

    % Use the parameters found via minimization to calculate the 
    % probability distribution from the Ising model.
    % The binomial prefactor is explained in Macke 2011. 
    Q_unnormalized = binom.*exp(states*param_list);

    % Normalize this Q.
    % This gives us a probability distribution Q(k), the probability of
    % finding a population with count k.
    Q_Ising = Q_unnormalized/sum(Q_unnormalized);

    % Function to raise the probability dist to the power 1/T.
    % Non-trivial.
    P_beta = prob_dist_power(P,1/T,binom);
    Q_beta = prob_dist_power(Q,1/T,binom);
    PEIF_beta = prob_dist_power(P_EIF,1/T,binom);
    Ising_beta = prob_dist_power(Q_Ising,1/T,binom);

    % Find non-zeros (because we will take the log of P_beta).
    ind = find(P_beta~=0);

    % Calculate the heat capacity. HC = var(log(p(x)))/n. 
    c(j) = sum(P_beta(ind).*(log2(P_beta(ind))-log2(binom(ind))).^2) -...
    sum(P_beta(ind).*(log2(P_beta(ind))-log2(binom(ind))))^2;
    c(j) = c(j)/N(j);
    
    ind = find(Q_beta~=0);

    clf(j) = sum(Q_beta(ind).*(log2(Q_beta(ind))-log2(binom(ind))).^2) -...
    sum(Q_beta(ind).*(log2(Q_beta(ind))-log2(binom(ind))))^2;
    clf(j) = clf(j)/N(j);
    
    ind = find(PEIF_beta~=0);

    ceif(j) = sum(PEIF_beta(ind).*(log2(PEIF_beta(ind))-log2(binom(ind))).^2) -...
    sum(PEIF_beta(ind).*(log2(PEIF_beta(ind))-log2(binom(ind))))^2;
    ceif(j) = ceif(j)/N(j);
    
    ind = find(Ising_beta~=0);

    cising(j) = sum(Ising_beta(ind).*(log2(Ising_beta(ind))-log2(binom(ind))).^2) -...
    sum(Ising_beta(ind).*(log2(Ising_beta(ind))-log2(binom(ind))))^2;
    cising(j) = cising(j)/N(j);
end

% Plot the line T=1.
% plot(N,c,'b')
% hold on
% plot(N,clf,'r')
% plot(N,ceif,'b--')
% plot(N,cising,'k')

%%

heat_cap = [[c';ceif';clf';cising'],[N';N';N';N'],[ones(size(c'));2*ones(size(c'));3*ones(size(c'));4*ones(size(c'))]];