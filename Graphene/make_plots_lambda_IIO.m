% plot
clear all

SavePlots=1;
mode = 4; % choose functions
model_name = 'BFPV'; %or ALMA
% model_name = 'ALMA';

if mode == 1
    name = 'expcos';
    load('data_IIO_expcos_eps20_VACVAC_35to40_gbar0025_ALMA.mat');
%     load('data_IIO_expcos_eps20_VACVAC_35to40_gbar01_ALMA.mat');
%     load('data_IIO_expcos_eps20_VACVAC_35to40_gbar1_ALMA.mat');
%     load('data_IIO_expcos_eps20_VACVAC_32to37_gbar0025_BFPV.mat');
%     load('data_IIO_expcos_eps20_VACVAC_32to37_gbar01_BFPV.mat');
%     load('data_IIO_expcos_eps20_VACVAC_32to37_gbar1_BFPV.mat');
end
if mode == 2
    name = 'cos2';
%     load('data_IIO_cos2_eps20_VACVAC_35to40_gbar0025_ALMA.mat');
%     load('data_IIO_cos2_eps20_VACVAC_35to40_gbar01_ALMA.mat');
%     load('data_IIO_cos2_eps20_VACVAC_35to40_gbar1_ALMA.mat');
%     load('data_IIO_cos2_eps20_VACVAC_32to37_gbar0025_BFPV.mat');
%     load('data_IIO_cos2_eps20_VACVAC_32to37_gbar01_BFPV.mat');
    load('data_IIO_cos2_eps20_VACVAC_32to37_gbar1_BFPV.mat');
end
if mode == 4
    name = 'cos4';
%     load('data_IIO_cos4_eps20_VACVAC_35to40_gbar0025_ALMA.mat');
%     load('data_IIO_cos4_eps20_VACVAC_35to40_gbar01_ALMA.mat');
%     load('data_IIO_cos4_eps20_VACVAC_35to40_gbar1_ALMA.mat'); 
%     load('data_IIO_cos4_eps20_VACVAC_35to45_gbar1_ALMA.mat');
%     load('data_IIO_cos4_eps20_VACVAC_32to37_gbar0025_BFPV.mat');
%     load('data_IIO_cos4_eps20_VACVAC_32to37_gbar01_BFPV.mat');
%     load('data_IIO_cos4_eps20_VACVAC_32to37_gbar1_BFPV.mat');
    load('data_IIO_cos4_eps20_VACVAC_34to40_gbar1_BFPV.mat');
end
if mode == 8
end


fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('lambda_st = %g  lambda_end = %g\n',lambda(1),lambda(end));
fprintf('N_theta = %d N = %d N_lambda = %d N_eps = %d\n',N_theta,N,N_lambda,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

Iu_norm1 = Iu_norm(end,:); Iw_norm1 = Iw_norm(end,:);
BU_norm1 = BU_norm(end,:); BU_norm11 = BU_norm(1,:);
Qu_norm1 = Qu_norm(end,:); Sw_norm1 = Sw_norm(end,:);
Qu_norm0 = Qu_norm(1,:); Sw_norm0 = Sw_norm(1,:);
BU_ratio = BU_norm1./BU_norm11;


% figure(1);
% norm_max = max([max(Iu_norm1),max(Iw_norm1)]);
% plot(lambda,Iu_norm1,'b-o',lambda,Iw_norm1,'g-*');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|I^{(u)}|_2$ and $|I^{(w)}|_2$','interpreter','latex');
% title('$|I^{(u)}|_2$ and $|I^{(w)}|_2$ versus $\lambda$','interpreter','latex');
% ll = legend('$|I^{(u)}|_2$','$|I^{(w)}|_2$');
% set(ll,'FontSize',16,'interpreter','latex');
% 
% figure(2);
% norm_max = max([max(Qu_norm1),max(Sw_norm1)]);
% plot(lambda,Qu_norm1,'b-o',lambda,Sw_norm1,'g-*');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|Q|_2$ and $|S|_2$','interpreter','latex');
% title('$|Q|_2$ and $|S|_2$ versus $\lambda$','interpreter','latex');
% ll = legend('$|Q|_2$','$|S|_2$');
% set(ll,'FontSize',16,'interpreter','latex');
% 
% figure(3);
% plot(lambda,BU_ratio,'b-o');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('ratio');
% title('$|BU|_2$ ratio versus $\lambda$','interpreter','latex');

figure(4);
norm_max = max([max(Qu_norm1),max(Qu_norm0)]);
plot(lambda,Qu_norm1,'b-o',lambda,Qu_norm0,'r-*');
xlabel('$\lambda$','interpreter','latex','FontSize',20);
ylabel('$|Q|_2$','interpreter','latex','FontSize',20);
title('$|Q|_2$ versus $\lambda$','interpreter','latex','FontSize',20);
ll = legend('$|Q|_2$ at $\varepsilon=\bar{g}/5$','$|Q|_2$ at $\varepsilon=0$');
set(ll,'FontSize',20,'interpreter','latex','Location','best');

figure(5);
norm_max = max([max(Sw_norm1),max(Sw_norm0)]);
plot(lambda,Sw_norm1,'b-o',lambda,Sw_norm0,'r-*');
xlabel('$\lambda$','interpreter','latex','FontSize',20);
ylabel('$|S|_2$','interpreter','latex','FontSize',20);
title('$|S|_2$ versus $\lambda$','interpreter','latex','FontSize',20);
ll = legend('$|S|_2$ at $\varepsilon=\bar{g}/5$','$|S|_2$ at $\varepsilon=0$');
set(ll,'FontSize',20,'interpreter','latex','Location','best');


if(SavePlots==1)
%     filename = sprintf('fig_index_%s_eps%.0f%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(2,filename,'epsc');
    filename = sprintf('fig_LGSPR_IIO_Qu_shift_%s_gbar%.0f_%s_%s',name,a*1000,IN,model_name);
    saveas(4,filename,'epsc');
    filename = sprintf('fig_LGSPR_IIO_Sw_shift_%s_gbar%.0f_%s_%s',name,a*1000,IN,model_name);
    saveas(5,filename,'epsc');
end