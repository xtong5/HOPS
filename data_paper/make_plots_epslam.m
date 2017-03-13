% plot
clear all
close all
SavePlots = 1;

mode = 1; % choose functions

if mode == 1
    name = 'expcos';
     load('data_expcos5_eps_VACAg.mat');
%     load('data_expcos5_eps_WATERAg.mat');
%     load('data_expcos10_eps_VACAg.mat');
%     load('data_expcos10_eps_WATERAg.mat');
%     load('data_expcos100_eps_VACAg.mat');
%     load('data_expcos100_eps_WATERAg.mat');
end
if mode == 2
    name = 'cos2';
    load('data_cos25_eps_VACAg.mat');
%     load('data_cos25_eps_WATERAg.mat');
%     load('data_cos210_eps_VACAg.mat');
%     load('data_cos210_eps_WATERAg.mat');
%     load('data_cos2100_eps_VACAg.mat');
%     load('data_cos2100_eps_WATERAg.mat');
end
if mode == 4
    name = 'cos4';
%     load('data_cos45_eps_VACAg.mat');
%     load('data_cos45_eps_WATERAg.mat');
%     load('data_cos410_eps_VACAg.mat');
    load('data_cos410_eps_WATERAg.mat');
%     load('data_cos4100_eps_VACAg.mat');
%     load('data_cos4100_eps_WATERAg.mat');
end

fprintf('-------------\n');
fprintf('a = %g  Eps_max = %g\n',a,Eps_max);
fprintf('N_theta = %d N = %d N_lamb = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

figure(1);
contourf(lambda,epsvec,Gn_U_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|\widetilde{U}|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

figure(2);
contourf(lambda,epsvec,Gn_W_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|\widetilde{W}|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

if(SavePlots==1)
    filename = sprintf('fig_U_%s_epslam%s%s%.0f',name,IN,OUT,a/Eps_max);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_W_%s_epslam%s%s%.0f',name,IN,OUT,a/Eps_max);
    saveas(2,filename,'epsc');
end