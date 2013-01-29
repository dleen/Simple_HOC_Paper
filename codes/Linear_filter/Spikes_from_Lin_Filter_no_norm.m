function [mu,rho,P_Th,cv,double_count,spike_times,n] =...
 Spikes_from_Lin_Filter_no_norm(r_est,bin_size,N,gen_method,dt,num_of_bins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate spikes from firing rate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise spike train matrix
n = zeros(num_of_bins,N);
coeff_of_var = zeros(1,N);

%%%%%%%%%%%%%%%%%%%%
% Coin flip method %
%%%%%%%%%%%%%%%%%%%%
if strcmp(gen_method,'coin flip')
    temp=r_est*dt;
    %t=zeros(length(r_est),N);
    for j=1:N
        l=2;
        for k=1:num_of_bins
            for i=round((k-1)*(bin_size/dt)+1):round(k*(bin_size/dt))
                if(temp(i)>rand)
                    n(k,j)=n(k,j)+1;
                    t(l,j)=i*dt;
                    l=l+1;
                end
            end
        end
        max_len(j)=length(find(t(:,j)>0));
        t_temp=t(1:max_len(j),j);
        coeff_of_var(j) = std(diff(t_temp))/mean(diff(t_temp));
    end
    spike_times=t(1:min(max_len),:);
    double_count = 100*length(find(n>1))/(N*length(n));
    max(n)
    %n = histc(spike_times,0:bin_size:min(max(spike_times)));
    % Replace multiple spikes in a bin with a single spike (if necessary
    n(n>1) = 1;
    Th_spike_stack = n;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Maths for neuroscientists algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(gen_method,'maths for neuro')
    %t=zeros(length(r_est),N);
    max_len=zeros(N,1);
    temp=r_est*dt;
    r_len=length(r_est);
    for j=1:N
        i=2; k=2; 
        while(i<r_len)
            y0=exprnd(1);
            i=i-1; 
            temp_sum=0;
            while(temp_sum<y0 && i<r_len)
                temp_sum=temp_sum+temp(i);
                i=i+1;
            end
            %if(((i-1)*dt - t(k-1,j))>3)
            t(k,j)=(i-1)*dt;
            k=k+1;
            %end
            if(k>r_len), break; end
        end        
        max_len(j)=length(find(t(:,j)>0));
        t_temp=t(1:max_len(j),j);
        coeff_of_var(j) = std(diff(t_temp))/mean(diff(t_temp));
    end
    %spike_times=t(1:max(max_len),:);
    spike_times=t(1:min(max_len),:);
    n = histc(spike_times,0:bin_size:min(max(spike_times)));
    max(n)
    %n_binned = histc(spike_times,0:bin_size:max(max(spike_times)));
    % Check for multiple spikes in a bin
    double_count = 100*length(find(n>1))/(N*length(n));
    % Replace multiple spikes in a bin with a single spike (if necessary
    %max(n)
    % Call the final result Th_spike_stack
    Th_spike_stack = n;
    %Th_spike_stack(Th_spike_stack>1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Math for neuro from the book %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(gen_method,'mfn')
    n_ind = length(r_est);
    %t=zeros(length(r_est),N);
    max_len=zeros(N,1);
    y_vdt = r_est*dt;
    for j=1:N
        ind = 1;
        thres = exprnd(1);
        y_int = 0;
        k=1;
        while ( ind < n_ind )
            while ( (y_int < thres) && (ind < n_ind) )
                y_int = y_int + y_vdt(ind);
                ind = ind + 1;
            end;
            if ( y_int >= thres )
                t(k,j)=(ind-1)*dt;
                k=k+1;
                
                thres = exprnd(1);
                y_int = 0;
            end;
        end;
        max_len(j)=length(find(t(:,j)>0));
        t_temp=t(1:max_len(j),j);
        coeff_of_var(j) = std(diff(t_temp))/mean(diff(t_temp));
    end
    spike_times=t(1:min(max_len),:);
    n = histc(spike_times,0:bin_size:min(max(spike_times)));
    %n_binned = histc(spike_times,0:bin_size:max(max(spike_times)));
    % Check for multiple spikes in a bin
    double_count = 100*length(find(n>1))/(N*length(n));
    % Replace multiple spikes in a bin with a single spike (if necessary
    max(n)
    % Call the final result Th_spike_stack
    Th_spike_stack = n;
    %Th_spike_stack(Th_spike_stack>1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version of Poisson spike generator found on internet %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(gen_method,'peche')
    %t=zeros(length(r_est),N);
    for j=1:N
        [T,tvec]=PechePourPoisson(r_est,dt);
        max_len(j)=length(T);
        t(1:max_len(j),j)=T;
         coeff_of_var(j) = std(diff(T))/...
            mean(diff(T));
    end
    spike_times = t(1:min(max_len),:);
    n = histc(spike_times,0:bin_size:max(max(spike_times)));
        double_count = 100*length(find(n>1))/(N*length(n));
    % Replace multiple spikes in a bin with a single spike (if necessary)
    n(n>1) = 1;
    Th_spike_stack = n;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usual calculations for mean, correlations etc. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temp vector to hold the number of occurences.
occurr = zeros(N+1,1);
% Sum across the rows to get the population spike count.
holder = sum(Th_spike_stack,2); 

% In this loop we count the number of times each population spike count
% occurs. E.g. If holder = (1,1,1,4,3,3,7,2)^t.
% ind = find(holder == 3) will give the indices where the number 3 can 
% be found: ind = (5,6) and by taking the length of ind we get the 
% number of times 3 occurs. length(ind) = 2. This is the number of
% occurrences of the spike count 3.
for i=0:N
    tempind = find(holder == i); % ^^^^
    occurr(i+1) = length(tempind); % count the number of occurrences
end

% Finally, return the population spike count probability. This P(k) gives
% the probability that a state with k spikes will occur.
P_Th = occurr/sum(occurr);

% Statistics of the DG.
mu = mean(mean(Th_spike_stack)); % mean firing rate.
sigma_mat = cov(Th_spike_stack); % find the covariance of the output.
%rho = mean(sigma_mat(1,2:N))/mean(diag(sigma_mat));
rho=((sum(sum(sigma_mat))-sum(diag(sigma_mat)))/(N*(N-1)))/...
    mean(diag(sigma_mat));
% ^^^ defined in the paper.
cv = mean(coeff_of_var);
end