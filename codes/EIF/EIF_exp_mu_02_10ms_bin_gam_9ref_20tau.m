%clear all; clf;

if(exist('./EIF_statistics_1ms.mexmaci64','file')~=3)
    mex ./mex_source/mex_main.cpp...
        ./mex_source/LIF_spike.cpp...
        ./mex_source/LIF_gen_spike_matrix.cpp...
        ./mex_source/DG_gen_spike_matrix.cpp...
        ./mex_source/EIF_gen_spike_matrix.cpp...
        ./mex_source/Macke_figures.cpp...
        -lm -lgsl -lgslcblas -output EIF_statistics_10ms_ref_10
end

%Number of neurons
N = 100;

%Number of trials
NT = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.5, rho = 0.05 %
%%%%%%%%%%%%%%%%%%%%%%%%
vals_05 = zeros(4,NT);
P_EIF_05 = zeros(N+1,NT);

tic
%matlabpool(4)
%parfor i=1:NT
for i=1:NT
    [vals_05(:,i),P_EIF_05(:,i)]=EIF_statistics_1ms(-12,0.33,6.23,N);
end
toc

P_err_05 = std(P_EIF_05,0,2)/sqrt(NT);
vals_err_05 = std(vals_05,0,2)/sqrt(NT);

P_mean_05 = mean(P_EIF_05,2);
vals_mean_05 = mean(vals_05,2);

%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.5, rho = 0.1 %
%%%%%%%%%%%%%%%%%%%%%%%
vals_1 = zeros(4,NT);
P_EIF_1 = zeros(N+1,NT);


%parfor i=1:NT
for i=1:NT
    [vals_1(:,i),P_EIF_1(:,i)]=EIF_statistics_1ms(-12,0.54,6.23,N);
end

P_err_1 = std(P_EIF_1,0,2)/sqrt(NT);
vals_err_1 = std(vals_1,0,2)/sqrt(NT);

P_mean_1 = mean(P_EIF_1,2);
vals_mean_1 = mean(vals_1,2);

%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.5, rho = 0.25%
%%%%%%%%%%%%%%%%%%%%%%%
vals_25 = zeros(4,NT);
P_EIF_25 = zeros(N+1,NT);


%parfor i=1:NT
for i=1:NT
    [vals_25(:,i),P_EIF_25(:,i)]=EIF_statistics_1ms(-12,0.842,6.23,N);
end
%matlabpool close

P_err_25 = std(P_EIF_25,0,2)/sqrt(NT);
vals_err_25 = std(vals_25,0,2)/sqrt(NT);

P_mean_25 = mean(P_EIF_25,2);
vals_mean_25 = mean(vals_25,2);