% plot
clear all

mode = 1000; % choose functions

if mode == 1
    load('data_expcos_eps.mat');
%     load('data_expcos10.mat');
%     load('data_expcos100.mat');
%     load('data_expcos5.mat');
end
if mode == 2
    load('data_cos252.mat');
%     load('data_cos210.mat');
%     load('data_cos2100.mat');
%     load('data_cos25.mat');
end
if mode == 4
    load('data_cos410.mat');
end
if mode == 8
    load('data_cos8.mat');
end

load('data_NoPer.mat');

fprintf('-------------\n');
fprintf(' a = %g  b = %g\n',a,b);
fprintf('N_theta = %d N = %d M = %d\n',N_theta,N,M);
fprintf('\n');

% U_norm1 = U_norm(:,2);
% W_norm1 = W_norm(:,2);
% BU_norm1 = BU_norm(:,2);


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