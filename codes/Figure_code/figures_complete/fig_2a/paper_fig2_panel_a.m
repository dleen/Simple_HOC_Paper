clear all;
%clf

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 1;

warning('off')

% Number of neurons used in simulation
N_range = [8,32,64,100];

% P = numerical
P_list = [];
% Q = Ising
Q_list = [];
N_list = [];
N_size = [];
% Just a list of ones to keep 
% track of which model is which when
% we use this data in another program
% EIF = 1
% Ising = 2
Model_EIF = [];
Model_Ising = [];

for N=N_range
    N_list = [N_list; (0:N)'];
    N_size = [N_size; N*ones(N+1,1)];
    
    Model_EIF = [Model_EIF;ones(N+1,1)];
    Model_Ising = [Model_Ising;2*ones(N+1,1)];
    
    % Import the data to be used.
    % Data file has the form:
    % p(x_0) p(x_1) ... p(x_(N-1)) mu rho
    M{1} = importdata(strcat('../numerical_data/fig2mu/fig_2mu_',...
        int2str(N),'_',num2str(0.1),'.dat'),' ');
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
    % bigger explanation needed
    param_list_init = [0;0];


    % Take the mean of the input p(x).
    % i.e we average over the trials
    % p(x_1) p(x_2) ... p(x_N) mu rho % trial number 1
    % p(x_1) p(x_2) ... p(x_N) mu rho % trial number 2
    P = mean(PP,1)';

    % plotting thing
    P_list = [P_list; P];
    
    % Calculate the means and moments i.e. the quantities
    % sum(k*P_k) and sum(k^2*P_k).
    mean_feature = P'*states;

    % Choose the algorithm to use.
    if alg_choice == 1
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
    
    Q_list = [Q_list; Q];
end

% Easy format for plotting in R.
fig_2a_R =...
 [[P_list;Q_list],[N_size;N_size],[N_list;N_list],[Model_EIF;Model_Ising]];