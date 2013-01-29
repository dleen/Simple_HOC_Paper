cols = winter(3);

% p1=plot(0:N,P_DG_05,'color',cols(1,:));
% hold on
% p2=plot(0:N,P_DG_1,'color',cols(2,:));
% p3=plot(0:N,P_DG_25,'color',cols(3,:));

p1=semilogy(0:N,P_LF_mean_05,'color',cols(1,:));
hold on
p2=semilogy(0:N,P_LF_mean_1,'color',cols(2,:));
p3=semilogy(0:N,P_LF_mean_25,'color',cols(3,:));

%axis('tight')
axis([0 100 1e-4 1])
l1=legend([p1 p2 p3],'\rho = 0.05','\rho = 0.1','\rho = 0.25');
set(l1,'fontsize',16)

xlabel('Pop spike count','fontsize',16)
ylabel('P(n)','fontsize',16)

%title('Dichotomized Gaussian \mu=0.5','fontsize',18)