clear all;

trials=80;
j=1;

matlabpool(4)

for N=2:2:100
    tic
    N
    mu   = zeros(trials,1);
    rho  = zeros(trials,1);
    P_LF = zeros(trials,N+1);
    
    parfor i=1:trials
        [mu(i),rho(i),P_LF(i,:),~,~]=...
            EIF_filter(-61.35,0.565,6.19,N,'mfn','l',10,'none');
    end
    
    M{j} = {mu,rho,P_LF};
    j=j+1;
    toc
end

matlabpool close