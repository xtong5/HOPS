% plot
clear all
close all

mode = 1; % choose functions

if mode == 1
%     load('data_expcos5_eps.mat');
%     load('data_expcos10_eps.mat');
    load('data_expcos100_eps.mat');
end
if mode == 2
    load('data_cos25_eps.mat');
%     load('data_cos210_eps.mat');
%     load('data_cos2100_eps.mat');
end
if mode == 4
    load('data_cos410_eps.mat');
%     load('data_cos4100_eps.mat');    
end
if mode == 8
    load('data_cos810_eps.mat');
%     load('data_cos810_eps.mat');
end


fprintf('-------------\n');
fprintf('Eps = %g  a = %g  b = %g\n',epsvec(end),a,b);
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,size(epsvec,2));
fprintf('\n');


figure(1);
contourf(lambda,epsvec,U_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|U|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

figure(2);
contourf(lambda,epsvec,W_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|W|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

figure(3);
contourf(lambda,epsvec,BU_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|BU|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');
