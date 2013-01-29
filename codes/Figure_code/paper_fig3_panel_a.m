clear all;
%clf

DKL = zeros(49,1);
DJS = DKL;

%rho_range = [0.05,0.1,0.25];
rho_range = [0.05];

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 1;

j=1;
k=1;

warning('off')

%M = importdata('./Lin_filt_data_01_01_lin_unnormed.mat');

for N=4:2:100
    %for k=1:length(rho_range)
        % Import the data to be used.
        % Data file has the form:
        % p(x_0) p(x_1) ... p(x_(N-1)) mu rho
        MM{1} = importdata(strcat('./numerical_data/fig2mu/fig_2mu_',...
            int2str(N),'_',num2str(rho_range(k)),'.dat'),' ');
        % PP contains the probability distributions.
        % Each row represents a single run at mean mu and corr rho.
        PP_EIF = MM{1}(:,3:N+3);
        mu_EIF = MM{1}(:,N+4);
        rho_EIF = MM{1}(:,N+5);
        
        %data=mean(cell2mat(M{j+1}));
        %PP_LF = data(1,3:N+3);
        %mu_LF = data(1,1);
        %rho_LF= data(1,2);
        
        P = mean(PP_EIF,1)';
        %Q = PP_LF';
        
        Q_th = linear_theor_func(-60,0.184,6.19,N,10);
        %Q_th = linear_theor_func(-61.35,0.565,6.19,N,10);
        Q_th = Q_th';
        
        
        binom = zeros(N+1,1);
        for i=0:N
            binom(i+1) = nchoosek(N,i);
        end

        % Find non-zero entries as we will be taking the log.
        ind_p = find(P ~= 0);
        %ind_q = find(Q ~= 0);
        ind_q_th = find(Q_th ~= 0);

        % Calculate the DKL.
        %DKL(j,k) = -Q(ind_q)'*(log2(Q(ind_q))-log2(binom(ind_q))) +...
        %            P(ind_p)'*(log2(P(ind_p))-log2(binom(ind_p)));
        % That is only the DKL for maxent models.
        
        %DKL(j,k) = P(ind_q)'*log(P(ind_q)./Q(ind_q));
                
        %Md = (1/2)*(P+Q);
        
        %DJS(j,k) = (1/2)*(P(ind_p)'*log(P(ind_p)./Md(ind_p))+...
        %    Q(ind_q)'*log(Q(ind_q)./Md(ind_q)));
        
        %DJSl(j,k) = DJS(j,k)/log(N);
        
        
        Md_th = (1/2)*(P+Q_th);
        
        DJS_th(j,k) = (1/2)*(P(ind_p)'*log(P(ind_p)./Md_th(ind_p))+...
            Q_th(ind_q_th)'*log(Q_th(ind_q_th)./Md_th(ind_q_th)));
        
        DJSl_th(j,k) = DJS_th(j,k)/log(N);
        
    %end
    j=j+1;
end

%pkl3=plot(2:2:100,DKL,'b')
%hold on
%%
%pjs3=plot(4:2:100,DJSl,'b')
hold on
pjs4=plot(4:2:100,DJSl_th,'r')