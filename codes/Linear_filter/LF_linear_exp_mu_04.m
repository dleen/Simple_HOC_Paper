% DG values with standard errors

%Number of neurons
N = 100;

%Number of trials
NT = 3;

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.2, rho = 0.05 %
%%%%%%%%%%%%%%%%%%%%%%%%
mu_LF_05 = zeros(NT,1);
rho_LF_05 = zeros(NT,1);
P_LF_05 = zeros(N+1,NT);

for i=1:NT
    [mu_LF_05(i),rho_LF_05(i),P_LF_05(:,i),cv,dbl]=...
        EIF_filter(-57.99,0.18,6.295,N,'mfn','l',10);
end

P_LF_err_05 = std(P_LF_05,0,2)/sqrt(NT);
m_LF_err_05 = std(mu_LF_05)/sqrt(NT);
r_LF_err_05 = std(rho_LF_05)/sqrt(NT);

P_LF_mean_05 = mean(P_LF_05,2);
m_LF_mean_05 = mean(mu_LF_05);
r_LF_mean_05 = mean(rho_LF_05);

%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.4, rho = 0.1 %
%%%%%%%%%%%%%%%%%%%%%%%
mu_LF_1 = zeros(NT,1);
rho_LF_1 = zeros(NT,1);
P_LF_1 = zeros(N+1,NT);

for i=1:NT
    [mu_LF_1(i),rho_LF_1(i),P_LF_1(:,i),cv,dbl]=...
        EIF_filter(-54.5,0.45,6.295,N,'mfn','l',10,'normed');
end

P_LF_err_1 = std(P_LF_1,0,2)/sqrt(NT);
m_LF_err_1 = std(mu_LF_1)/sqrt(NT);
r_LF_err_1 = std(rho_LF_1)/sqrt(NT);

P_LF_mean_1 = mean(P_LF_1,2);
m_LF_mean_1 = mean(mu_LF_1);
r_LF_mean_1 = mean(rho_LF_1);

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.2, rho = 0.25 %
%%%%%%%%%%%%%%%%%%%%%%%%
mu_LF_25 = zeros(NT,1);
rho_LF_25 = zeros(NT,1);
P_LF_25 = zeros(N+1,NT);

for i=1:NT
    [mu_LF_25(i),rho_LF_25(i),P_LF_25(:,i),cv,dbl]=...
        EIF_filter(-60.13,0.735,6.295,N,'mfn','l',10);
end

P_LF_err_25 = std(P_LF_25,0,2)/sqrt(NT);
m_LF_err_25 = std(mu_LF_25)/sqrt(NT);
r_LF_err_25 = std(rho_LF_25)/sqrt(NT);

P_LF_mean_25 = mean(P_LF_25,2);
m_LF_mean_25 = mean(mu_LF_25);
r_LF_mean_25 = mean(rho_LF_25);