% plot
clear all

SavePlots=0;
mode = 1; % choose functions

if mode == 1  
%     load('data_expcos_eps1_VACAg.mat');   
%     load('data_expcos_eps10_VACAg.mat');
    load('data_expcos_eps10_VACAg_test1.mat');   
%     load('data_expcos_eps20_VACAg.mat');   
%     load('data_expcos_eps40_VACAg.mat');   
%     load('data_expcos_eps1_WATERAg.mat');   
%     load('data_expcos_eps10_WATERAg.mat');   
%     load('data_expcos_eps20_WATERAg.mat');   
%     load('data_expcos_eps40_WATERAg.mat');    
end
if mode == 2
%     load('data_cos2_eps1_VACAg.mat');   
%     load('data_cos2_eps10_VACAg.mat');   
%     load('data_cos2_eps20_VACAg.mat');  
    load('data_cos2_eps20_VACAg_test.mat');   
%     load('data_cos2_eps40_VACAg.mat');   
%     load('data_cos2_eps1_WATERAg.mat');   
%     load('data_cos2_eps10_WATERAg.mat');   
%     load('data_cos2_eps20_WATERAg.mat');   
%     load('data_cos2_eps40_WATERAg.mat');    
end
if mode == 4
    load('data_cos4_eps1_VACAg.mat');   
%     load('data_cos4_eps10_VACAg.mat');   
%     load('data_cos4_eps20_VACAg.mat');   
%     load('data_cos4_eps1_WATERAg.mat');   
%     load('data_cos4_eps10_WATERAg.mat');   
%     load('data_cos4_eps20_WATERAg.mat');   
end
if mode == 8
    load('data_cos8_eps1_VACAg.mat');   
%     load('data_cos8_eps10_VACAg.mat');   
%     load('data_cos8_eps20_VACAg.mat');   
%     load('data_cos8_eps1_WATERAg.mat');   
%     load('data_cos8_eps10_WATERAg.mat');   
%     load('data_cos8_eps20_WATERAg.mat');   
end


fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

U_norm1 = U_norm(end,:);
W_norm1 = W_norm(end,:);
U_norm1_dno = U_norm_dno(end,:);
W_norm1_dno = W_norm_dno(end,:);
% Qu_norm1 = Qu_norm(end,:);
% Sw_norm1 = Sw_norm(end,:);
% BU_norm1 = BU_norm(end,:);
% BU_norm11 = BU_norm(1,:);
% BU_ratio = BU_norm1./BU_norm11;
Gn_U_norm1 = Gn_U_norm(end,:);
Gn_W_norm1 = Gn_W_norm(end,:);
Gn_U_norm1_dno = Gn_U_norm_dno(end,:);
Gn_W_norm1_dno = Gn_W_norm_dno(end,:);

lambda_crit_plot = lambda_crit*ones(M,1)';
% norm_max = max([max(U_norm1),max(W_norm1),max(BU_norm1)]);
% norm_max = max([max(U_norm1),max(W_norm1),max(Qu_norm1),max(Sw_norm1),max(BU_norm1)]);

% norm_max = max([max(Qu_norm1),max(Sw_norm1)]);
% yy = linspace(0,norm_max,M)';
% 
% figure(1);
% % plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',...
% %     lambda,Gn_U_norm1,'c-d',lambda,Gn_W_norm1,'y-+',lambda_crit_plot,yy,'r--');
% % plot(lambda,Qu_norm1,'c-d',lambda,Sw_norm1,'y-+',lambda,BU_norm1,'b-o',lambda_crit_plot,yy,'r--');
% % plot(lambda,U_norm1,'c-d',lambda,W_norm1,'y-+',lambda,BU_norm1,'b-o',lambda_crit_plot,yy,'r--');
% plot(lambda,U_norm1,'c-d',lambda,W_norm1,'y-+',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|Qu|_2$ and $|Sw|_2$','interpreter','latex');
% title('$|Qu|_2$ and $|Sw|_2$ versus $\lambda$','interpreter','latex');
% % legend('|Qu|_2','|Sw|_2','|BU|_2','lambda_c');
% legend('|Qu|_2','|Sw|_2','lambda_c');

norm_max = max(U_norm1);
yy = linspace(0,norm_max,M)';
figure(111);
plot(lambda,U_norm1,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$','interpreter','latex');
title('$|U|_2$ versus $\lambda$','interpreter','latex');
legend('|U|_2','lambda_c');

% norm_max = max(Qu_norm1);
% yy = linspace(0,norm_max,M)';
% figure(112);
% plot(lambda,Qu_norm1,'c-d',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|Qu|_2$','interpreter','latex');
% title('$|Qu|_2$ versus $\lambda$','interpreter','latex');
% legend('|Qu|_2','lambda_c');

norm_max = max(Gn_U_norm1);
yy = linspace(0,norm_max,M)';
figure(113);
plot(lambda,Gn_U_norm1,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|GU|_2$','interpreter','latex');
title('$|GU|_2$ versus $\lambda$','interpreter','latex');
legend('|GU|_2','lambda_c');

norm_max = max(W_norm1);
yy = linspace(0,norm_max,M)';
figure(121);
plot(lambda,W_norm1,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|W|_2$','interpreter','latex');
title('$|W|_2$ versus $\lambda$','interpreter','latex');
legend('|W|_2','lambda_c');

% norm_max = max(Sw_norm1);
% yy = linspace(0,norm_max,M)';
% figure(122);
% plot(lambda,Sw_norm1,'c-d',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|Sw|_2$','interpreter','latex');
% title('$|Sw|_2$ versus $\lambda$','interpreter','latex');
% legend('|Sw|_2','lambda_c');

norm_max = max(Gn_W_norm1);
yy = linspace(0,norm_max,M)';
figure(123);
plot(lambda,Gn_W_norm1,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|GW|_2$','interpreter','latex');
title('$|GW|_2$ versus $\lambda$','interpreter','latex');
legend('|GW|_2','lambda_c');

% norm_max = max(BU_norm1);
% yy = linspace(0,norm_max,M)';
% figure(13);
% plot(lambda,BU_norm1,'c-d',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|BU|_2$','interpreter','latex');
% title('$|BU|_2$ versus $\lambda$','interpreter','latex');
% legend('|BU|_2','lambda_c');


% figure(2);
% plot(lambda,BU_ratio,'b-o');
% xlabel('$\lambda$','interpreter','latex');
% % ylabel('$|Qu|_2$ and $|Sw|_2$','interpreter','latex');
% ylabel('ratio');
% title('$|BU|_2$ ratio versus $\lambda$','interpreter','latex');

% eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
% yy= linspace(0,eps_max,N_eps)';

% figure(3);
% plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
%     lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$\epsilon$','interpreter','latex');
% title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
% legend('epsilon_u','-Re[epsilon_w]','lambda_c');


norm_max = max([max(Gn_U_norm1),max(Gn_W_norm1)]);
yy = linspace(0,norm_max,M)';
figure(3);
hh = gca;
plot(lambda,Gn_U_norm1,'b-o',lambda,Gn_W_norm1,'g-*',lambda_crit_plot,yy,'r--');

xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$','interpreter','latex');
title('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
set(hh,'FontSize',16);
ll = legend('$|\widetilde{U}|_2$','$|\widetilde{W}|_2$','$\lambda_c$');
set(ll,'FontSize',16,'interpreter','latex');

norm_max = max([max(U_norm1),max(W_norm1)]);
yy = linspace(0,norm_max,M)';
figure(4);
hh = gca;
plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$ and $|W|_2$','interpreter','latex');
title('$|U|_2$ and $|W|_2$ versus $\lambda$','interpreter','latex');
set(hh,'FontSize',16);
ll = legend('$|U|_2$','$|W|_2$','$\lambda_c$');
set(ll,'FontSize',16,'interpreter','latex');


norm_max = max(U_norm1_dno);
yy = linspace(0,norm_max,M)';
figure(51);
plot(lambda,U_norm1_dno,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$','interpreter','latex');
title('$|U|_2$ DNO versus $\lambda$','interpreter','latex');
legend('|U|_2','lambda_c');

norm_max = max(Gn_U_norm1_dno);
yy = linspace(0,norm_max,M)';
figure(52);
plot(lambda,Gn_U_norm1_dno,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|GU|_2$','interpreter','latex');
title('$|GU|_2$ DNO versus $\lambda$','interpreter','latex');
legend('|GU|_2','lambda_c');

norm_max = max(W_norm1_dno);
yy = linspace(0,norm_max,M)';
figure(61);
plot(lambda,W_norm1_dno,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|W|_2$','interpreter','latex');
title('$|W|_2$ DNO versus $\lambda$','interpreter','latex');
legend('|W|_2','lambda_c');

norm_max = max(Gn_W_norm1_dno);
yy = linspace(0,norm_max,M)';
figure(62);
plot(lambda,Gn_W_norm1_dno,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|GW|_2$','interpreter','latex');
title('$|GW|_2$ DNO versus $\lambda$','interpreter','latex');
legend('|GW|_2','lambda_c');


if(SavePlots==1)
    filename = sprintf('fig_UWlam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_BUratio_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(2,filename,'epsc');
    filename = sprintf('fig_index_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(3,filename,'epsc');
end