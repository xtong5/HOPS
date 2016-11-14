% plot
clear all

mode = 4; % choose functions

load('data_NoPer.mat');
U_norm_NoPer = U_norm;
W_norm_NoPer = W_norm;
BU_norm_NoPer = BU_norm;

if mode == 1
    load('data_expcos100.mat');
    U_norm_100 = U_norm;
    W_norm_100 = W_norm;
    BU_norm_100 = BU_norm;
    load('data_expcos10.mat');
    U_norm_10 = U_norm;
    W_norm_10 = W_norm;
    BU_norm_10 = BU_norm;
end
if mode == 2
    load('data_cos2100.mat');
    U_norm_100 = U_norm;
    W_norm_100 = W_norm;
    BU_norm_100 = BU_norm;
    load('data_cos210.mat');
    U_norm_10 = U_norm;
    W_norm_10 = W_norm;
    BU_norm_10 = BU_norm;
end
if mode == 4
    load('data_cos4100.mat');
    U_norm_100 = U_norm;
    W_norm_100 = W_norm;
    BU_norm_100 = BU_norm;
    load('data_cos410.mat');
    U_norm_10 = U_norm;
    W_norm_10 = W_norm;
    BU_norm_10 = BU_norm;
end
if mode == 8
    load('data_cos8100.mat');
    U_norm_100 = U_norm;
    W_norm_100 = W_norm;
    BU_norm_100 = BU_norm;
    load('data_cos810.mat');
    U_norm_10 = U_norm;
    W_norm_10 = W_norm;
    BU_norm_10 = BU_norm;
end
if mode == 0
%     load('data_square.mat');
%     U_norm_100 = U_norm;
%     W_norm_100 = W_norm;
%     BU_norm_100 = BU_norm;
    load('data_squareEps1.mat');
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
end


fprintf('-------------\n');
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('N_theta = %d N = %d M = %d\n',N_theta,N,M);
fprintf('\n');

if mode >0

lambda_crit_plot = lambda_crit*ones(M,1);
norm_max = max([max(U_norm_NoPer),max(U_norm_100),max(U_norm_10)]);
yy = linspace(0,norm_max,M)';

figure(1);
plot(lambda,U_norm_NoPer,'b-o',lambda,U_norm_100,'g-*',...
    lambda,U_norm_10,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$','interpreter','latex');
title('$|U|_2$ versus $\lambda$','interpreter','latex');
legend('|No|_2','|100|_2','|10|_2','lambda_c');

norm_max = max([max(W_norm_NoPer),max(W_norm_100),max(W_norm_10)]);
yy = linspace(0,norm_max,M)';

figure(2);
plot(lambda,W_norm_NoPer,'b-o',lambda,W_norm_100,'g-*',...
    lambda,W_norm_10,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|W|_2$','interpreter','latex');
title('$|W|_2$ versus $\lambda$','interpreter','latex');
legend('|No|_2','|100|_2','|10|_2','lambda_c');

norm_max = max([max(BU_norm_NoPer),max(BU_norm_100),max(BU_norm_10)]);
yy = linspace(0,norm_max,M)';

figure(3);
plot(lambda,BU_norm_NoPer,'b-o',lambda,BU_norm_100,'g-*',...
    lambda,BU_norm_10,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|BU|_2$','interpreter','latex');
title('$|BU|_2$ versus $\lambda$','interpreter','latex');
legend('|No|_2','|100|_2','|10|_2','lambda_c');
end
    

% eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
% yy= linspace(0,eps_max,M)'; 
% 
% figure(2);
% plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
%     lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$\epsilon$','interpreter','latex');
% title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
% legend('epsilon_u','-Re[epsilon_w]','lambda_c');