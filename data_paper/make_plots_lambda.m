% plot
clear all
close all

SavePlots = 0;

mode = 4; % choose functions

if mode == 1
    name = 'expcos';
%     load('data_expcos5_eps_VACAg.mat');
%     load('data_expcos5_eps_WATERAg.mat');
%     load('data_expcos10_eps_VACAg.mat');
%     load('data_expcos10_eps_WATERAg.mat');
%     load('data_expcos100_eps_VACAg.mat');
%     load('data_expcos100_eps_WATERAg.mat');
    load('data_expcos15_eps_VACAg.mat');
%     load('data_expcos15_eps_WATERAg.mat');
end
if mode == 2
    name = 'cos2';
    load('data_cos25_eps_VACAg.mat');
%     load('data_cos25_eps_WATERAg.mat');
%     load('data_cos210_eps_VACAg.mat');
%     load('data_cos210_eps_WATERAg.mat');
%     load('data_cos2100_eps_VACAg.mat');
%     load('data_cos2100_eps_WATERAg.mat');
%     load('data_cos215_eps_VACAg.mat');
    load('data_cos215_eps_WATERAg.mat');
end
if mode == 4
    name = 'cos4';
%     load('data_cos45_eps_VACAg.mat');
%     load('data_cos45_eps_WATERAg.mat');
%     load('data_cos410_eps_VACAg.mat');
%     load('data_cos410_eps_WATERAg.mat');
%     load('data_cos4100_eps_VACAg.mat');
%     load('data_cos4100_eps_WATERAg.mat');
%     load('data_cos45_eps_VACAgN20.mat');
%     load('data_cos45_eps_VACAgN24.mat');
%     load('data_cos45_eps_VACAgNt96N16.mat');
%     load('data_cos45_eps_VACAgNt96N20.mat');
%     load('data_cos45_eps_VACAgNt96N24.mat');
%     load('data_cos45_eps_WATERAgN20.mat');
%     load('data_cos45_eps_WATERAgN24.mat');
%     load('data_cos45_eps_WATERAgNt96N16.mat');
%     load('data_cos45_eps_WATERAgNt96N20.mat');
%     load('data_cos45_eps_WATERAgNt96N24.mat');
%     load('data_cos415_eps_VACAg.mat');
    load('data_cos415_eps_WATERAg.mat');
end


fprintf('-------------\n');
fprintf('a = %g  Eps_max = %g\n',a,Eps_max);
fprintf('N_theta = %d N = %d N_lamb = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

U_norm1 = U_norm(end,:);
W_norm1 = W_norm(end,:);
Gn_U_norm1 = Gn_U_norm(end,:);
Gn_W_norm1 = Gn_W_norm(end,:);

lambda_crit_plot = lambda_crit*ones(M,1)';
% norm_max = max([max(U_norm1),max(W_norm1));
norm_max = max([max(Gn_U_norm1),max(Gn_W_norm1)]);

yy = linspace(0,norm_max,M)';

figure(1);
% plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',...
%     lambda,Gn_U_norm1,'c-d',lambda,Gn_W_norm1,'y-+',lambda_crit_plot,yy,'r--');
% plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',lambda_crit_plot,yy,'r--');
plot(lambda,Gn_U_norm1,'b-o',lambda,Gn_W_norm1,'g-*',lambda_crit_plot,yy,'r--');

xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$','interpreter','latex');
title('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
legend('|U|_2','|W|_2','lambda_c');

eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
yy= linspace(0,eps_max,M)';

figure(2);
plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
legend('epsilon_u','-Re[epsilon_w]','lambda_c');

if(SavePlots==1)
    filename = sprintf('fig_UWlam_%s_%s%s%.0f',name,IN,OUT,...
        Eps_max/a*100);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_index_%s_eps%s%s%.0f',name,IN,OUT,...
        Eps_max/a*100);
    saveas(2,filename,'epsc');
end
