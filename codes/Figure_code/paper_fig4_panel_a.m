clear all;

rho_range = [0.05,0.1,0.25];

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 1;

j=1;

warning('off')

N=100;

P = [];

for k=1:3
    % Import the data to be used.
    % Data file has the form:
    % p(x_0) p(x_1) ... p(x_(N-1)) mu rho
    M{1} = importdata(strcat('./numerical_data/fig2mu/fig_2mu_',...
        int2str(N),'_',num2str(rho_range(k)),'.dat'),' ');
    % PP contains the probability distributions.
    % Each row represents a single run at mean mu and corr rho.
    PP = M{1}(:,3:N+3);
    mu = M{1}(:,N+4);
    rho = M{1}(:,N+5);

    P = [P;mean(PP,1)'];
end

[~,~,~,P_DG_05]=DG_statistics(0.132,-1.28,100);
[~,~,~,P_DG_1]=DG_statistics(0.242,-1.28,100);
[~,~,~,P_DG_25]=DG_statistics(0.502,-1.28,100);

Q = [P_DG_05;P_DG_1;P_DG_25];

rho_vals = [0.05*ones(101,1);0.1*ones(101,1);0.25*ones(101,1)];

nums = [(0:100)';(0:100)';(0:100)'];

EIF_DG_comp = [[P;Q],[nums;nums],[rho_vals;rho_vals],[ones(size(P));2*ones(size(P))]];