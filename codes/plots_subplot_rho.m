figure(1)
subplot(1,3,1)
title('Correlation coeff    \rho = 0.05','fontsize',15)
p11=semilogy(0:N,P_DG_mean_05,'color','r');
hold on
p12=semilogy(0:N,P_LF_mean_05,'color','b');
xlabel('Population spike count')
ylabel('P(n)','fontsize',16)

subplot(1,3,2)
title('Correlation coeff    \rho = 0.1','fontsize',15)
p21=semilogy(0:N,P_DG_mean_1,'color','r');
hold on
p22=semilogy(0:N,P_LF_mean_1,'color','b');
xlabel('Population spike count')

subplot(1,3,3)
title('Correlation coeff    \rho = 0.25','fontsize',15)
p31=semilogy(0:N,P_DG_mean_25,'color','r');
hold on
p32=semilogy(0:N,P_LF_mean_25,'color','b');
xlabel('Population spike count')

l1=legend([p11 p12],'DG','Lin Filter');
l2=legend([p21 p22],'DG','Lin Filter');
l3=legend([p31 p32],'DG','Lin Filter');


set(l1,'fontsize',14)
set(l2,'fontsize',14)
set(l3,'fontsize',14)

