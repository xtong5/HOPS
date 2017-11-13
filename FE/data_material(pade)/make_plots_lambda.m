% plot
clear all

mode = 8; % choose functions

if mode == 1
    load('data_expcos10_eps_WATERAg.mat');
    %load('data_cos410_eps_WATERAg.mat');
end
if mode == 2
    load('data_cos45_eps_WATERAg.mat');
end
if mode == 4
    load('data_cos410.mat');
end
if mode == 8
    load('data_cos810_eps_WATERAg.mat');
end

% load('data_NoPer.mat');

fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

U_norm1 = U_norm(end,:);
%U_norm1 = U_norm(1,:);
W_norm1 = W_norm(end,:);
BU_norm1 = BU_norm(end,:);


lambda_crit_plot = lambda_crit*ones(M,1);
norm_max = max([max(U_norm1),max(W_norm1),max(BU_norm1)]);
yy = linspace(0,norm_max,M)';

figure(1);
plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',...
    lambda,BU_norm1,'c-d',lambda_crit_plot,yy,'r--');
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