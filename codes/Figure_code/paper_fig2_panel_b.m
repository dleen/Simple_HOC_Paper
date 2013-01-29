clear all;
clf

DKL = zeros(49,3);
DJS = DKL;

rho_range = [0.05,0.1,0.25];

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 1;

j=1;

warning('off')

for N=4:2:100
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

        % Find non-zero entries as we will be taking the log.
        ind_p = find(P ~= 0);
        ind_q = find(Q ~= 0);

        % Calculate the DKL.
        %DKL(j,k) = -Q(ind_q)'*(log2(Q(ind_q))-log2(binom(ind_q))) +...
        %            P(ind_p)'*(log2(P(ind_p))-log2(binom(ind_p)));
        %DKL(j,k) = DKL(j,k)/log(N);
        
        Md_th = (1/2)*(P+Q);
        
        DJS(j,k) = (1/2)*(P(ind_p)'*log(P(ind_p)./Md_th(ind_p))+...
            Q(ind_q)'*log(Q(ind_q)./Md_th(ind_q)));
        
        DJS(j,k) = DJS(j,k)/log(N);        
    end
    j=j+1;
end