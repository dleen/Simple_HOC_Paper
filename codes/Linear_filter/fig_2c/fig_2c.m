clear all;

trials=10;
j=1;

for N=2:2:100
    mu   = zeros(trials,1);
    rho  = zeros(trials,1);
    P_LF = zeros(trials,N+1);
    
    for i=1:trials
        [mu(i),rho(i),P_LF(i,:),~,~]=...
            EIF_filter(-56.15,0.183,6.295,N,'mfn','l',10,'normed');
    end
    
    M{j} = {mu,rho,P_LF};
    j=j+1;
end