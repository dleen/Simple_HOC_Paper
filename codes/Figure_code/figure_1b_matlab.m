clear all;
tic
clf
CM = hsv(4);

% Choose between Matlab optimization = 1
% or minimize.m minimizaçtion function = 2.
alg_choice = 1;

% Number of neurons.
N = 3;
% Number of repeated trials at fixed mu, rho.
avnum = 5;

% Load data from a file. 
% The file.dat has the form:
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 1
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 2
% We keep mu and rho fixed over trials. We have multiple
% trials so that we can get an average and smooth out
% some of the inherent noise.
M{1} = importdata(strcat('./numerical_data/fig2mu/fig_2mu_',...
    int2str(N),'_0.05','.dat'),' ');

% Generate a list of binomial coefficients.
binom = zeros(N+1,1);
for i=0:N
    binom(i+1) = nchoosek(N,i);
end

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

% An initial guess:
param_list_init = [0;0];
param_list = param_list_init;

%DKL = zeros(length(index),4);
%DKL_sorted = DKL;

%    DKL_averaged = zeros(floor(length(index)/avnum),4);
%   rho_averaged = zeros(floor(length(index)/avnum),4);
%   DKL_error = rho_averaged;

for j=1:4
    % PP is a matrix where each row is a prob dist for a given
    % value of mu and rho.
    [B,MIX] = sort(M{j}(:,1));

    PP  = M{j}(MIX,3:N+3);
    mu  = M{j}(MIX,N+4);
    rho = M{j}(MIX,N+5);
    
    index = find(~isnan(rho));
    rho = rho(index);
    mu = mu(index);
    PP = PP(index,:);

    DKL = zeros(length(index));
    
    for k=1:length(index)
        P = (PP(k,:))';
        % Calculate the means and moments i.e. the quantities
        % sum(k*P_k) and sum(k^2*P_k).
        mean_feature = P'*states;

        if alg_choice == 1
            % MATLAB Optimization Toolbox minimization function.
            % Uses Optimization toolbox unconstrained minimization function.
            options = optimset('GradObj','on','LargeScale','on',...
            'Display','off',...
            'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
            % We set GradObj = on as we supply the gradient.
            % We set LargeScale = on as that is the algorithm that uses the
            % user supplied gradient.
            % Display just shows some accuracy output.
            % The final options are tolerances etc.

            param_list_init = param_list;
            % The minimization function
            [param_list,fval,output] =...
                fminunc(@(x)neg_log_like_binom...
                (x,states,mean_feature,P,binom),param_list_init,options);
        elseif alg_choice == 2
            % Eric's minimization function, based on conjugate gradient.
            param_list_init = param_list;
            param_list = minimize(param_list_init,...
                'neg_log_like_binom',200,...
                states,mean_feature,P,binom);
        end

        % Use the parameters found via minimization to calculate the
        % probability distribution from the Ising model.
        Q_unnormalized = binom.*exp(states*param_list);

        % Normalize this Q.
        % This gives us a probability distribution Q(k), the probability of
        % finding a population with count k.
        Q = Q_unnormalized/sum(Q_unnormalized);

        % Find non-zero entries as we will be taking the log.
        ind_p = find(P ~= 0);
        ind_q = find(Q ~= 0);
        %ind_p = find(P > 1e-14);
        %ind_q = find(Q > 1e-14);

        % Calculate the DKL.
        %DKL(k,j) = -Q(ind_q)'*log2(Q(ind_q)) + P(ind_p)'*log2(P(ind_p));
        %DKL(k) = -Q(ind_q)'*log2(Q(ind_q)) + P(ind_p)'*log2(P(ind_p));
        DKL(k) = -Q(ind_q)'*(log2(Q(ind_q))-log2(binom(ind_q))) +...
                    P(ind_p)'*(log2(P(ind_p))-log2(binom(ind_p)));

    end

    % Sort rho to get it in increasing order.
    [rho_sorted,IX] = sort(rho);
    % Match corresponding DKLs.
    %DKL_sorted(:,j) = DKL(IX,j);
    DKL = DKL(IX);


    % Average over avnum values to get a smoother plot.
%     for i=1:floor(length(index)/avnum)
%         DKL_averaged(i,j) = mean(DKL_sorted((avnum*(i-1)+1):avnum*i,j));
%         rho_averaged(i,j) = mean(rho_sorted((avnum*(i-1)+1):avnum*i));
%         DKL_error(i,j) = (1/sqrt(avnum))*std(DKL_sorted((avnum*(i-1)+1):(avnum*i),j));
%     end

    av_end = floor(length(index)/avnum);
    DKL_averaged = zeros(av_end,1);
    DKL_error = DKL_averaged;
    rho_averaged = DKL_averaged;
    for i=1:av_end
        DKL_averaged(i) = mean(DKL((avnum*(i-1)+1):avnum*i));
        rho_averaged(i) = mean(rho_sorted((avnum*(i-1)+1):avnum*i));
        DKL_error(i) = (1/sqrt(avnum))*std(DKL((avnum*(i-1)+1):(avnum*i)));
    end

    % Sort once more.
    [rho_final,IX] = sort(rho_averaged);
    DKL_final = DKL_averaged(IX);
    error_final = DKL_error(IX);

    [h1(j),hp(j)]=boundedline(rho_final,DKL_final,error_final);
    %plot(rho_final,DKL_final,'color',CM(j,:))
    hold on
end


xlabel('Correlation coefficient','fontsize',16)
ylabel('KL-divergence \Delta_h','fontsize',16)

axis([0 0.3 0 10e-3])
for i=1:4
    set(h1(i),'color',CM(i,:));
    set(hp(i),'facecolor',CM(i,:))
    alpha(0.3)
    set(get(get(hp(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
end

leg1b = legend('\mu = 0.02','\mu = 0.05','\mu = 0.1','\mu = 0.2','location','NorthWest');
legend('boxoff')
set(leg1b,'fontsize',16)

title('EIF 10ms bin size','fontsize',18)
box on;

% set(hp2,'facecolor',CM(1,:))
% alpha(0.3)
% set(h2,'color',CM(1,:))
% %[h3,hp3]=boundedline(rho_final(:,3),DKL_final(:,3),error_final(:,3));
% 
% set(get(get(hp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
% 
% legend('\mu = 0.02','\mu = 0.05','\mu = 0.1','\mu = 0.2')

toc