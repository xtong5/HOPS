% plot
clear all
close all


load('data_NoPer.mat');
U_norm_NoPer = U_norm;
W_norm_NoPer = W_norm;
BU_norm_NoPer = BU_norm;

load('data_squareEps1pade.mat');
% load('data_squareEps1.mat');
U_norm_Eps1 = U_norm;
W_norm_Eps1 = W_norm;
BU_norm_Eps1 = BU_norm;
    
lambda_crit_plot = lambda_crit*ones(M,1);
norm_max = max([max(U_norm_NoPer),max(U_norm_Eps1)]);
yy = linspace(0,norm_max,M)';

figure(1);
plot(lambda,U_norm_NoPer,'b-o',lambda,U_norm_Eps1,'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$','interpreter','latex');
title('$|U|_2$ versus $\lambda$','interpreter','latex');
legend('|No|_2','|Eps1|_2','lambda_c');

norm_max = max([max(W_norm_NoPer),max(W_norm_Eps1)]);
yy = linspace(0,norm_max,M)';

figure(2);
plot(lambda,W_norm_NoPer,'b-o',lambda,W_norm_Eps1,'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|W|_2$','interpreter','latex');
title('$|W|_2$ versus $\lambda$','interpreter','latex');
legend('|No|_2','|Eps1|_2','lambda_c');

norm_max = max([max(BU_norm_NoPer),max(BU_norm_Eps1)]);
yy = linspace(0,norm_max,M)';

figure(3);
plot(lambda,BU_norm_NoPer,'b-o',lambda,BU_norm_Eps1,'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|BU|_2$','interpreter','latex');
title('$|BU|_2$ versus $\lambda$','interpreter','latex');
legend('|No|_2','|Eps1|_2','lambda_c');

