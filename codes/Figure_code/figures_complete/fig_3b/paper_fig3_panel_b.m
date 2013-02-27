clear all;
%clf

DKL = zeros(49,3);
DJS = DKL;

rho_range = [0.05,0.1,0.25];
rho_vals_DG = [0.132,0.242,0.502];
rho_vals_LF = [0.215,0.37,0.66];
mu_vals = [-60.31,-60.78,-63.6];

warning('off')

for k=1:length(rho_range)
    jj=1;
    
    lambda_DG = rho_vals_DG(k);

    lambda_LF = rho_vals_LF(k);

    gamma_dg = -1.28;


    %%%%%%%%%%%%%%
    % LIF values %
    %%%%%%%%%%%%%%
    gamma=mu_vals(k);
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

    r_est = r_nl_full;

    for i=1:(num_of_bins-floor(bin_size/dt))
        SS(i) =...
        sum(r_est(((i-1)*floor(bin_size/dt)+1):(i*floor(bin_size/dt))))*dt;
    end

    s_dg = -3:0.01:5;
    s_lf = s_dg;

    [plf,plf_s] = ksdensity(SS,s_lf);

    EIF = 1-exp(-plf_s);
    EIF = max(0,EIF);

    P_LF_func = @(kk,N) nchoosek(N,kk)*...
        sum(plf.*((1-EIF).^(N-kk)).*((EIF).^kk))*(plf_s(2)-plf_s(1));
    
    for N=4:2:100
        % Import the data to be used.
        % Data file has the form:
        % p(x_0) p(x_1) ... p(x_(N-1)) mu rho
        MM{1} = importdata(strcat('../numerical_data/fig2mu/fig_2mu_',...
            int2str(N),'_',num2str(rho_range(k)),'.dat'),' ');
        % PP contains the probability distributions.
        % Each row represents a single run at mean mu and corr rho.
        PP_EIF = MM{1}(:,3:N+3);
        mu_EIF = MM{1}(:,N+4);
        rho_EIF = MM{1}(:,N+5);
        
        P = mean(PP_EIF,1)';
        
        clear Q
        for i=0:N
            Q(i+1) = P_LF_func(i,N);
        end
        Q=Q';
        
        binom = zeros(N+1,1);
        for i=0:N
            binom(i+1) = nchoosek(N,i);
        end

        % Find non-zero entries as we will be taking the log.
        ind_p = find(P ~= 0);
        ind_q = find(Q ~= 0);
        
        Md_th = (1/2)*(P+Q);
        
        DJS(jj,k) = (1/2)*(P(ind_p)'*log(P(ind_p)./Md_th(ind_p))+...
            Q(ind_q)'*log(Q(ind_q)./Md_th(ind_q)));
        
        DJS(jj,k) = DJS(jj,k)/log(N);
        
        jj=jj+1;
        N
    end
end