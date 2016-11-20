% plot
clear all
close all

mode = 8; % choose functions

if mode == 1
%     load('data_expcos100_eps_VACAg.mat');
    load('data_expcos100_eps_VACAu.mat');
%     load('data_expcos100_eps_WATERAg.mat');
%     load('data_expcos100_eps_WATERAu.mat');
%     load('data_expcos10_eps_VACAg.mat');
%     load('data_expcos10_eps_VACAu.mat');
%     load('data_expcos10_eps_WATERAg.mat');
%     load('data_expcos10_eps_WATERAu.mat');
%     load('data_expcos5_eps_VACAg.mat');
%     load('data_expcos5_eps_VACAu.mat');
%     load('data_expcos5_eps_WATERAg.mat');
%     load('data_expcos5_eps_WATERAu.mat');
%     load('data_expcos52_eps_VACAg.mat');
%     load('data_expcos52_eps_VACAu.mat');
%     load('data_expcos52_eps_WATERAg.mat');
%     load('data_expcos52_eps_WATERAu.mat');   
end
if mode == 2
%     load('data_cos2100_eps_VACAg.mat');
%     load('data_cos2100_eps_VACAu.mat');
%     load('data_cos2100_eps_WATERAg.mat');
%     load('data_cos2100_eps_WATERAu.mat');
%     load('data_cos210_eps_VACAg.mat');
%     load('data_cos210_eps_VACAu.mat');
%     load('data_cos210_eps_WATERAg.mat');
%     load('data_cos210_eps_WATERAu.mat');
%     load('data_cos25_eps_VACAg.mat');
%     load('data_cos25_eps_VACAu.mat');
%     load('data_cos25_eps_WATERAg.mat');
%     load('data_cos25_eps_WATERAu.mat');
%     load('data_cos252_eps_VACAg.mat');
%     load('data_cos252_eps_VACAu.mat');
%     load('data_cos252_eps_WATERAg.mat');
    load('data_cos252_eps_WATERAu.mat');
end
if mode == 4
%     load('data_cos4100_eps_VACAg.mat');
%     load('data_cos4100_eps_VACAu.mat');
%     load('data_cos4100_eps_WATERAg.mat');
%     load('data_cos4100_eps_WATERAu.mat');
%     load('data_cos410_eps_VACAg.mat');
%     load('data_cos410_eps_VACAu.mat');
%     load('data_cos410_eps_WATERAg.mat');
%     load('data_cos410_eps_WATERAu.mat');
%     load('data_cos45_eps_VACAg.mat');
%     load('data_cos45_eps_VACAu.mat');
%     load('data_cos45_eps_WATERAg.mat');
%     load('data_cos45_eps_WATERAu.mat');
end
if mode == 8
%     load('data_cos8100_eps_VACAg.mat');
%     load('data_cos8100_eps_VACAu.mat');
%     load('data_cos8100_eps_WATERAg.mat');
%     load('data_cos8100_eps_WATERAu.mat');
%     load('data_cos810_eps_VACAg.mat');
%     load('data_cos810_eps_VACAu.mat');
%     load('data_cos810_eps_WATERAg.mat');
%     load('data_cos810_eps_WATERAu.mat');
%     load('data_cos85_eps_VACAg.mat');
%     load('data_cos85_eps_VACAu.mat');
%     load('data_cos85_eps_WATERAg.mat');
    load('data_cos85_eps_WATERAu.mat');
end


fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');


figure(1);
contourf(lambda,epsvec,U_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|U|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

figure(2);
contourf(lambda,epsvec,W_norm);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|W|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

figure(3);
contourf(lambda,epsvec,BU_norm);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|BU|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');
