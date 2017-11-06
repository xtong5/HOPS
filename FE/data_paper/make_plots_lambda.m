% plot
clear all
close all

SavePlots = 0;

mode = 4; % choose functions

if mode == 1
    name = 'expcos';
%     load('data_expcos_eps10_VACAg.mat');
%     load('data_expcos_eps10_WATERAg.mat');
%     load('data_expcos_eps20_VACAg.mat');
    load('data_expcos_eps20_WATERAg.mat');
end
if mode == 2
    name = 'cos2';
%     load('data_cos2_eps20_VACAg.mat');
%     load('data_cos2_eps20_WATERAg.mat');
%     load('data_cos2_eps10_VACAg.mat');
    load('data_cos2_eps10_WATERAg.mat');
end
if mode == 4
    name = 'cos4';
%     load('data_cos4_eps10_VACAg.mat');
%     load('data_cos4_eps10_WATERAg.mat');
%     load('data_cos4_eps20_VACAgNt128N16.mat');
    load('data_cos4_eps20_WATERAgNt128N16.mat');
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
hh = gca;
% plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',...
%     lambda,Gn_U_norm1,'c-d',lambda,Gn_W_norm1,'y-+',lambda_crit_plot,yy,'r--');
% plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',lambda_crit_plot,yy,'r--');
plot(lambda,Gn_U_norm1,'b-o',lambda,Gn_W_norm1,'g-*',lambda_crit_plot,yy,'r--');

xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$','interpreter','latex');
title('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
set(hh,'FontSize',16);
ll = legend('$|\widetilde{U}|_2$','$|\widetilde{W}|_2$','$\lambda_c$');
set(ll,'FontSize',16,'interpreter','latex');

eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
yy= linspace(0,eps_max,M)';

figure(2);
hh = gca;
plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$\epsilon_u$ and -Re[$\epsilon_w$] versus $\lambda$','interpreter','latex');
set(hh,'FontSize',16);
ll = legend('$\epsilon_u$','-Re[$\epsilon_w$]','$\lambda_c$');
set(ll,'FontSize',16,'interpreter','latex');

if(SavePlots==1)
    filename = sprintf('fig_UWlam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_index_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(2,filename,'epsc');
end
