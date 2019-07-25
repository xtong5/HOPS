% plot
clear all

SavePlots=0;
function_mode = 2; % choose functions

if function_mode == 1
    name = 'expcos';
%     load('data_expcos_eps10_VACAg.mat');
    load('data_expcos_eps20_VACAg.mat');
%     load('data_expcos_eps10_WATERAg.mat');
%     load('data_expcos_eps20_WATERAg.mat');
end
if function_mode == 2
    name = 'cos2';
%     load('data_cos2_eps10_VACAg.mat');
    load('data_cos2_eps20_VACAg.mat');
%     load('data_cos2_eps10_WATERAg.mat');
%     load('data_cos2_eps20_WATERAg.mat');   
end
if function_mode == 4
    name = 'cos4';
    load('data_cos4_eps10_VACAg.mat');   
%     load('data_cos4_eps20_VACAg_Nth128.mat');   
%     load('data_cos4_eps10_WATERAg.mat');   
%     load('data_cos4_eps20_WATERAg_Nth128.mat');   
end



fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

U_norm1 = U_norm(end,:);
W_norm1 = W_norm(end,:);
Qu_norm1 = Qu_norm(end,:);
Sw_norm1 = Sw_norm(end,:);
BU_norm1 = BU_norm(end,:);
BU_norm11 = BU_norm(1,:);
Gn_U_norm1 = Gn_U_norm(end,:);
Gn_W_norm1 = Gn_W_norm(end,:);
BU_ratio = BU_norm1./BU_norm11;



lambda_crit_plot = lambda_crit*ones(M,1)';

figure(1);
% norm_max = max([max(U_norm1),max(W_norm1),max(BU_norm1)]);
norm_max = max([max(U_norm1)]);
yy = linspace(0,norm_max,M)';
plot(lambda,U_norm1,'b-d');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$' ,'interpreter','latex');
title('$|U|_2$  versus $\lambda$','interpreter','latex');
legend('|U|_2');


% figure(2);
% norm_max = max([max(Gn_U_norm1),max(Gn_W_norm1)]);
% yy = linspace(0,norm_max,M)';
% plot(lambda,Gn_U_norm1,'c-d',lambda,Gn_W_norm1,'y-+',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$','interpreter','latex');
% title('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
% ll = legend('$|\widetilde{U}|_2$','$|\widetilde{W}|_2$','$\lambda_c$');
% set(ll,'FontSize',16,'interpreter','latex');
% 
% figure(3);
% norm_max = max([max(Qu_norm1),max(Sw_norm1)]);
% yy = linspace(0,norm_max,M)';
% plot(lambda,Qu_norm1,'c-d',lambda,Sw_norm1,'y-+',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|Qu|_2$ and $|Sw|_2$','interpreter','latex');
% title('$|Qu|_2$ and $|Sw|_2$ versus $\lambda$','interpreter','latex');
% legend('|Qu|_2','|Sw|_2','lambda_c');
% 
% figure(4);
% plot(lambda,BU_ratio,'b-o');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('ratio');
% title('$|BU|_2$ ratio versus $\lambda$','interpreter','latex');
% 
% % figure(5);
% % eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
% % yy= linspace(0,eps_max,M)';
% % plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
% %     lambda_crit_plot,yy,'r--');
% % xlabel('$\lambda$','interpreter','latex');
% % ylabel('$\epsilon$','interpreter','latex');
% % title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
% % legend('epsilon_u','-Re[epsilon_w]','lambda_c');

if(SavePlots==1)
    filename = sprintf('fig_IIO_UWlam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(2,filename,'epsc');
    filename = sprintf('fig_IIO_QuSwlam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(3,filename,'epsc');
    filename = sprintf('fig_IIO_BUratiolam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(4,filename,'epsc');
%     filename = sprintf('fig_index_%s_eps%.0f%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(5,filename,'epsc');
end

saveas(1,'LSPR_1','epsc');