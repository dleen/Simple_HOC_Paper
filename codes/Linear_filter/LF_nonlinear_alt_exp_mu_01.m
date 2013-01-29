% DG values with standard errors

%Number of neurons
N = 100;

%Number of trials
NT = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.1, rho = 0.05 %
%%%%%%%%%%%%%%%%%%%%%%%%
mu_LF_05 = zeros(NT,1);
rho_LF_05 = zeros(NT,1);
P_LF_05 = zeros(N+1,NT);

for i=1:NT
    [mu_LF_05(i),rho_LF_05(i),P_LF_05(:,i),cv,dbl]=...
    EIF_filter_nl_change(-60.6,0.31,5.15,N,'mfn','lnl_full',10,'normed');
    [mu_LF_05_none(i),rho_LF_05_none(i),P_LF_05_none(:,i),cv,dbl]=...
    EIF_filter_nl_change(-60.59,0.3,5.08,N,'mfn','lnl_full',10,'none');
end

P_LF_err_05 = std(P_LF_05,0,2)/sqrt(NT);
m_LF_err_05 = std(mu_LF_05)/sqrt(NT);
r_LF_err_05 = std(rho_LF_05)/sqrt(NT);

P_LF_mean_05 = mean(P_LF_05,2);
m_LF_mean_05 = mean(mu_LF_05);
r_LF_mean_05 = mean(rho_LF_05);

%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.1, rho = 0.1 %
%%%%%%%%%%%%%%%%%%%%%%%
mu_LF_1 = zeros(NT,1);
rho_LF_1 = zeros(NT,1);
P_LF_1 = zeros(N+1,NT);

for i=1:NT
    [mu_LF_1(i),rho_LF_1(i),P_LF_1(:,i),cv,dbl]=...
    EIF_filter_nl_change(-61.12,0.448,5.15,N,'mfn','lnl_full',10,'normed');
    [mu_LF_1(i),rho_LF_1(i),P_LF_1(:,i),cv,dbl]=...
    EIF_filter_nl_change(-61.11,0.425,5.08,N,'mfn','lnl_full',10,'none');    
end

P_LF_err_1 = std(P_LF_1,0,2)/sqrt(NT);
m_LF_err_1 = std(mu_LF_1)/sqrt(NT);
r_LF_err_1 = std(rho_LF_1)/sqrt(NT);

P_LF_mean_1 = mean(P_LF_1,2);
m_LF_mean_1 = mean(mu_LF_1);
r_LF_mean_1 = mean(rho_LF_1);

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.1, rho = 0.25 %
%%%%%%%%%%%%%%%%%%%%%%%%
mu_LF_25 = zeros(NT,1);
rho_LF_25 = zeros(NT,1);
P_LF_25 = zeros(N+1,NT);

for i=1:NT
    [mu_LF_25(i),rho_LF_25(i),P_LF_25(:,i),cv,dbl]=...
    EIF_filter_nl_change(-63.7,0.675,5.15,N,'mfn','lnl_full',10,'normed');
    [mu_LF_25_none(i),rho_LF_25_none(i),P_LF_25_none(:,i),cv,dbl]=...
    EIF_filter_nl_change(-63.96,0.657,5.08,N,'mfn','lnl_full',10,'none');
end

P_LF_err_25 = std(P_LF_25,0,2)/sqrt(NT);
m_LF_err_25 = std(mu_LF_25)/sqrt(NT);
r_LF_err_25 = std(rho_LF_25)/sqrt(NT);

P_LF_mean_25 = mean(P_LF_25,2);
m_LF_mean_25 = mean(mu_LF_25);
r_LF_mean_25 = mean(rho_LF_25);