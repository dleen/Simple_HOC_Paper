clear all;
clf

% Set the number of neurons used:
N = 96;

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 1;

% Import the data to be used.
% Data file has the form:
% p(x_0) p(x_1) ... p(x_(N-1)) mu rho
M{1} = importdata(strcat('./numerical_data/fig2mu/fig_2mu_',...
    int2str(N),'_0.05','.dat'),' ');

% PP contains the probability distributions.
% Each row represents a single run at mean mu and corr rho.
PP = M{1}(:,3:N+3);
mu = M{1}(:,N+4);
rho = M{1}(:,N+5);

% List the possible states, using the symmetry argument.
% Symmetry argument can be found in Macke 2011.
% Essentially says that the parameters in the Ising model h_i, J_ij
% using symmetry conditions can be reduced to just 2 parameters. The 
% symmetry condition means that we talk about population spike 
% counts rather than particular probabilities. i.e P(101101) becomes 
% just P(4) as 4 spikes occurred. 
% Here we list the possible states which are just 1...N and at the
% same time we list their squares as needed in the Ising model.
states = [(0:N)',(0:N)'.^2];

% Generate a list of binomial coefficients.
%warning('off');
binom = zeros(N+1,1);
for i=0:N
    binom(i+1) = nchoosek(N,i);
end

% An initial guess at the parameters to be minimized:
param_list_init = [0;0];


% Take the mean of the input p(x).
% i.e we average over the trials
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 1
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 2
P = mean(PP,1)';

% Calculate the means and moments i.e. the quantities
% sum(k*P_k) and sum(k^2*P_k).
mean_feature = P'*states;

% Choose the algorithm to use.
if alg_choice == 1
    % MATLAB Optimization Toolbox minimization function.
    % Uses Optimization toolbox unconstrained minimization function.
    options = optimset('GradObj','on','LargeScale','on',...
        'Display','final-detailed',...
        'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-15,'TolX',1e-14);
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
          (x,states,mean_feature,P,binom),param_list_init,options);
%     [param_list,fval,exitflag,output] =...
%           fsolve(@(x)neg_log_like_binom...
%           (x,states,mean_feature,P,binom),param_list_init,options);
elseif alg_choice == 2
    % "Minimize.m" function, based on conjugate gradient.
    param_list = minimize(param_list_init,...
          'neg_log_like_binom',200,...
          states,mean_feature,P,binom);
end  

% Use the parameters found via minimization to calculate the 
% probability distribution from the Ising model.
% The binomial prefactor is explained in Macke 2011. 
Q_unnormalized = binom.*exp(states*param_list);

% Normalize this Q.
% This gives us a probability distribution Q(k), the probability of
% finding a population with count k.
Q = Q_unnormalized/sum(Q_unnormalized);

% Find non-zero entries as we will be taking the log.
ind_p = find(P ~= 0);
ind_q = find(Q ~= 0);

% Calculate the DKL.
DKL = -Q(ind_q)'*(log2(Q(ind_q))-log2(binom(ind_q))) +...
            P(ind_p)'*(log2(P(ind_p))-log2(binom(ind_p)))


% Plot data.
if(N<50)
    bar([mean(PP)' Q])
elseif(N>=50 && N<200)
    p1=plot(0:N,P);
    hold on
    p2=plot(0:N,Q,'r');
    set(p1,'Color','blue','LineWidth',3);
    set(p2,'Color','red','LineWidth',3);
elseif(N==200)
    p1=plot(0:N,mean(PP));
    hold on
    p2=plot(0:N,Q,'r');
    set(p1,'Color','blue','LineWidth',3);
    set(p2,'Color','red','LineWidth',3);
end

% title('EIF Numerical Simulation','fontsize',18)
% fig2a_leg = legend('\rho = 0.05','\rho = 0.1','\rho = 0.25');
% set(fig2a_leg,'FontSize',16);
% xlabel('Population spike count i.e synchrony','fontsize',16)
% ylabel('Probability','fontsize',16)
% axis([0 N 1e-4 1]);
% hold off

% %Bar plot of p and p2 in layer 6    
set(gca,'FontSize',24)
title('P(n)')
xlabel('N=100 neurons')
axis([-Inf Inf 0 1])
%       bar_labels={'000','001','011','111'} ;
      bar_labels={} ;
    legend({'DG','Ising Fit'})
       set(gca, 'XTickLabel', bar_labels) 
axis([-0.1 Inf -Inf Inf])
