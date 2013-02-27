clear all;
close all;

warning('off')

N_range = [8,32,64,100];

P_list = [];
Q_list = [];
N_list = [];
N_size = [];
Model_EIF = [];
Model_Ising = [];

for N=N_range
    N_list = [N_list; (0:N)'];
    N_size = [N_size; N*ones(N+1,1)];
    
    Model_EIF = [Model_EIF;ones(N+1,1)];
    Model_Ising = [Model_Ising;2*ones(N+1,1)];


    lambda_DG = 0.242;

    lambda_LF = 0.363;

    gamma_dg = -1.28;


    %%%%%%%%%%%%%%
    % LIF values %
    %%%%%%%%%%%%%%
    gamma=-60.87;
    E0=gamma; % DC current
    tau_m=5; % membrane time constant
    tau_ref=3; % refractory period
    deltat=3; % EIF "AP sharpness"
    v_soft=-53; % EIF soft threshold
    v_reset=-60; % Reset value
    v_th=20; % EIF threshold
    sigma=4.87;
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

    pdg = normpdf(s_dg,0,sqrt(lambda_DG));
    pdg_s = s_dg; 

    % DG
    xx=(s_dg+gamma_dg)/sqrt(1-lambda_DG);
    erftest = erf(xx/sqrt(2));
    DG = 0.5*(1+erftest);

    EIF = 1-exp(-plf_s);
    EIF = max(0,EIF);


    P_DG_func = @(kk) nchoosek(N,kk)*...
        sum(pdg.*((1-DG).^(N-kk)).*((DG).^kk))*(pdg_s(2)-pdg_s(1));
    P_LF_func = @(kk) nchoosek(N,kk)*...
        sum(plf.*((1-EIF).^(N-kk)).*((EIF).^kk))*(plf_s(2)-plf_s(1));

    for i=0:N
        P(i+1) = P_DG_func(i);
        Q(i+1) = P_LF_func(i);
    end
   
    P_list = [P_list; P'];  
    Q_list = [Q_list; Q'];
end

fig_3a_R =...
 [[P_list;Q_list],[N_size;N_size],[N_list;N_list],[Model_EIF;Model_Ising]];