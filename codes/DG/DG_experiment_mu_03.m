% DG values with standard errors

%Number of neurons
N = 100;

%Number of trials
NT = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.3, rho = 0.05 %
%%%%%%%%%%%%%%%%%%%%%%%%
mu_05 = zeros(NT,1);
rho_05 = zeros(NT,1);
P_DG_05 = zeros(N+1,NT);

for i=1:NT
    [mu_05(i),rho_05(i),sigma,P_DG_05(:,i)]=DG_statistics(0.086,-0.525,N);
end

P_DG_err_05 = std(P_DG_05,0,2)/sqrt(NT);
m_DG_err_05 = std(mu_05)/sqrt(NT);
r_DG_err_05 = std(rho_05)/sqrt(NT);

P_DG_mean_05 = mean(P_DG_05,2);
m_DG_mean_05 = mean(mu_05);
r_DG_mean_05 = mean(rho_05);

%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.3, rho = 0.1 %
%%%%%%%%%%%%%%%%%%%%%%%
mu_1 = zeros(NT,1);
rho_1 = zeros(NT,1);
P_DG_1 = zeros(N+1,NT);

for i=1:NT
    [mu_1(i),rho_1(i),sigma,P_DG_1(:,i)]=DG_statistics(0.17,-0.525,N);
end

P_err_1 = std(P_DG_1,0,2)/sqrt(NT);
m_err_1 = std(mu_1)/sqrt(NT);
r_err_1 = std(rho_1)/sqrt(NT);

P_mean_1 = mean(P_DG_1,2);
m_mean_1 = mean(mu_1);
r_mean_1 = mean(rho_1);

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.3, rho = 0.25 %
%%%%%%%%%%%%%%%%%%%%%%%%
mu_25 = zeros(NT,1);
rho_25 = zeros(NT,1);
P_DG_25 = zeros(N+1,NT);

for i=1:NT
    [mu_25(i),rho_25(i),sigma,P_DG_25(:,i)]=DG_statistics(0.404,-0.525,N);
end

P_DG_err_25 = std(P_DG_25,0,2)/sqrt(NT);
m_DG_err_25 = std(mu_25)/sqrt(NT);
r_DG_err_25 = std(rho_25)/sqrt(NT);

P_DG_mean_25 = mean(P_DG_25,2);
m_DG_mean_25 = mean(mu_25);
r_DG_mean_25 = mean(rho_25);