% plot
clear all

SavePlots=0;
mode = 1; % choose functions
model_name = 'BFPV'; %or ALMA

if mode == 1
    name = 'expcos';
%     load('data_expcos_eps20_VACVAC_35to40_gbar1.mat');
%     load('data_expcos_eps20_VACVAC_35to40_gbar05.mat');
%     load('data_expcos_eps20_VACVAC_35to40_gbar005.mat');
%     load('data_expcos_eps20_VACVAC_35to40_gbar0025.mat');
%     load('data_expcos_eps20_VACVAC_32to37_gbar0025_BFPV.mat');
%     load('data_expcos_eps20_VACVAC_32to37_gbar005_BFPV.mat');
%     load('data_expcos_eps20_VACVAC_32to37_gbar05_BFPV.mat');
%     load('data_expcos_eps20_VACVAC_32to37_gbar1_BFPV.mat');
    load('data_expcos_eps20_VACVAC_34to40_gbar2_BFPV.mat');
%     load('data_expcos_eps20_VACVAC_32to37_gbar2_BFPV.mat');
%     load('data_expcos_eps20_VACVAC_30to40_gbar5_BFPV.mat');
%     load('data_expcos_eps20_VACVAC_25to40_gbar5_BFPV.mat');
end
if mode == 2
    name = 'cos2';
%     load('data_cos2_eps20_VACVAC_35to40_gbar1.mat');
%     load('data_cos2_eps20_VACVAC_35to40_gbar05.mat');
%     load('data_cos2_eps20_VACVAC_35to40_gbar005.mat');
%     load('data_cos2_eps20_VACVAC_35to40_gbar0025.mat');
%     load('data_cos2_eps20_VACVAC_32to37_gbar0025_BFPV.mat');
%     load('data_cos2_eps20_VACVAC_32to37_gbar005_BFPV.mat');
%     load('data_cos2_eps20_VACVAC_32to37_gbar05_BFPV.mat');
    load('data_cos2_eps20_VACVAC_32to37_gbar1_BFPV.mat');
end
if mode == 4
    name = 'cos4';
%     load('data_cos4_eps20_VACVAC_35to45_gbar1.mat');
    load('data_cos4_eps20_VACVAC_35to40_gbar1.mat');
%     load('data_cos4_eps20_VACVAC_35to40_gbar05.mat');
%     load('data_cos4_eps20_VACVAC_35to40_gbar005.mat');
%     load('data_cos4_eps20_VACVAC_35to40_gbar0025.mat');
%     load('data_cos4_eps20_VACVAC_32to37_gbar0025_BFPV.mat');
%     load('data_cos4_eps20_VACVAC_32to37_gbar005_BFPV.mat');
%     load('data_cos4_eps20_VACVAC_32to37_gbar05_BFPV.mat');
%     load('data_cos4_eps20_VACVAC_32to37_gbar1_BFPV.mat');
    load('data_cos4_eps20_VACVAC_34to40_gbar1_BFPV.mat');
end
if mode == 8
end


fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('lambda_st = %g  lambda_end = %g\n',lambda(1),lambda(end));
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

U_norm1 = U_norm(end,:);
W_norm1 = W_norm(end,:);
U_norm0 = U_norm(1,:);
W_norm0 = W_norm(1,:);
% BU_norm1 = BU_norm(end,:);
% BU_norm11 = BU_norm(1,:);
Gn_U_norm1 = Gn_U_norm(end,:);
Gn_W_norm1 = Gn_W_norm(end,:);
Gn_U_norm0 = Gn_U_norm(1,:);
Gn_W_norm0 = Gn_W_norm(1,:);
% BU_ratio = BU_norm1./BU_norm11;


figure(1);
norm_max = max([max(U_norm1),max(W_norm1)]);
plot(lambda,U_norm1,'b-o',lambda,W_norm1,'g-*');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$ and $|W|_2$','interpreter','latex');
title('$|U|_2$ and $|W|_2$ versus $\lambda$','interpreter','latex');
ll = legend('$|U|_2$','$|W|_2$');
set(ll,'FontSize',16,'interpreter','latex');

figure(2);
norm_max = max([max(Gn_U_norm1),max(Gn_W_norm1)]);
plot(lambda,Gn_U_norm1,'b-o',lambda,Gn_W_norm1,'g-*');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$','interpreter','latex');
title('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
ll = legend('$|\widetilde{U}|_2$','$|\widetilde{W}|_2$');
set(ll,'FontSize',16,'interpreter','latex');

% figure(3);
% plot(lambda,BU_ratio,'b-o');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('ratio');
% title('$|BU|_2$ ratio versus $\lambda$','interpreter','latex');


figure(4);
norm_max = max([max(Gn_U_norm1),max(Gn_U_norm0)]);
plot(lambda,Gn_U_norm1,'b-o',lambda,Gn_U_norm0,'r-*');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{U}|_2$','interpreter','latex');
title('$|\widetilde{U}|_2$ versus $\lambda$','interpreter','latex');
ll = legend('$|\widetilde{U}|_2$ at Eps=max','$|\widetilde{U}|_2$ at Eps=0');
set(ll,'FontSize',16,'interpreter','latex','Location','best');

figure(5);
norm_max = max([max(Gn_W_norm1),max(Gn_W_norm0)]);
plot(lambda,Gn_W_norm1,'b-o',lambda,Gn_W_norm0,'r-*');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{W}|_2$','interpreter','latex');
title('$|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
ll = legend('$|\widetilde{W}|_2$ at Eps=max','$|\widetilde{W}|_2$ at Eps=0');
set(ll,'FontSize',16,'interpreter','latex','Location','best');



if(SavePlots==1)
%     filename = sprintf('fig_index_%s_eps%.0f%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(2,filename,'epsc');
    filename = sprintf('fig_LGSPR_GU_shift_%s_gbar%.0f_%s_%s',name,a*1000,IN,model_name);
    saveas(4,filename,'epsc');
    filename = sprintf('fig_LGSPR_GW_shift_%s_gbar%.0f_%s_%s',name,a*1000,IN,model_name);
    saveas(5,filename,'epsc');
end