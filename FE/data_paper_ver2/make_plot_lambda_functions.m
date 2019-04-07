% plot
clear all
close all

SavePlots = 0;

mode = 4; % choose functions

load('data_expcos_eps10_VACAg.mat');U_expcos_eps10 = U_norm(end,:);
load('data_expcos_eps20_VACAg.mat');U_expcos_eps20 = U_norm(end,:);
load('data_cos2_eps10_VACAg.mat');U_cos2_eps10 = U_norm(end,:);
load('data_cos2_eps20_VACAg.mat');U_cos2_eps20 = U_norm(end,:);
load('data_cos4_eps10_VACAg.mat');U_cos4_eps10 = U_norm(end,:);
load('data_cos4_eps20_VACAgNt128N16.mat');U_cos4_eps20 = U_norm(end,:);


lambda_crit_plot = lambda_crit*ones(M,1)';
norm_max = max([max(U_expcos_eps10),max(U_expcos_eps20),...
    max(U_cos2_eps10),max(U_cos2_eps20),max(U_cos4_eps10),max(U_cos4_eps20)]);

yy = linspace(0,norm_max,M)';

figure(1);
hh = gca;
plot(lambda,U_expcos_eps10,'b-',lambda,U_expcos_eps20,'g-',lambda,U_cos2_eps10,'b-',lambda,U_cos2_eps20,'g-',...
    lambda,U_cos4_eps10,'b-',lambda,U_cos4_eps20,'g-');
hold on 

xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$','interpreter','latex');
title('Fields versus incident wavelength');
set(hh,'FontSize',16);
ll = legend('$|\widetilde{U}|_2$','$|\widetilde{W}|_2$','$\lambda_c$');
set(ll,'FontSize',16,'interpreter','latex');



