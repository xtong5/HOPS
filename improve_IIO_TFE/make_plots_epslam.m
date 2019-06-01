% plot
clear all
close all

mode = 1; % choose functions

if mode == 1  
    load('data_expcos_eps20_VACAg.mat');   
end
if mode == 2
end
if mode == 4
end
if mode == 8
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
