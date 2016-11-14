% plot
clear all
close all

mode = 2; % choose functions

if mode == 1
%     load('data_expcos100padeN12.mat');
    load('data_expcos10padeN12.mat');
end
if mode == 2
%     load('data_cos2100padeN12.mat');
    load('data_cos210padeN12.mat');
end
if mode == 4
%     load('data_cos4100padeN12.mat');
    load('data_cos410padeN12.mat');
end
if mode == 8
%     load('data_cos8100padeN12.mat');
    load('data_cos810padeN12.mat');
end

fprintf('-------------\n');
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('N_theta = %d N = %d M = %d\n',N_theta,N,M);
fprintf('\n');


lambda_crit_plot = lambda_crit*ones(M,1);
norm_max = max([max(U_norm),max(W_norm),max(BU_norm)]);
yy = linspace(0,norm_max,M)';

figure(1);
plot(lambda,U_norm,'b-o',lambda,W_norm,'g-*',...
    lambda,BU_norm,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$ and $|W|_2$','interpreter','latex');
title('$|U|_2$ and $|W|_2$ versus $\lambda$','interpreter','latex');
legend('|U|_2','|W|_2','|BU|_2','lambda_c');

eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
yy= linspace(0,eps_max,M)';

figure(2);
plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
legend('epsilon_u','-Re[epsilon_w]','lambda_c');