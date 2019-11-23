% PLOT_DELTA.m 
%
% Script to plot the determinant of DNO_TM
%
% XT 9/19

close all

SavePlots = 0;
L = 2*pi;
N_theta = 64;
M = 1001;
P = 4;

lambda_low = 340;
lambda_high = 350;
% lambda = linspace(0.3,0.8,M)';
lambda = linspace(lambda_low,lambda_high,M)';

k_0 = L./lambda;
n_u = zeros(N_lambda,1); epsilon_u = zeros(N_lambda,1);
n_w = zeros(N_lambda,1); epsilon_w = zeros(N_lambda,1);
IN = 'SILVER';
OUT = 'VACUUM';

for m = 1:N_lambda
    [n_u(m),epsilon_u(m)] = ri_perm(lambda(m),OUT);
    [n_w(m),epsilon_w(m)] = ri_perm(lambda(m),IN);
end
k_u = n_u.*k_0; 
k_w = n_w.*k_0;

%TM
tau2 = epsilon_u./epsilon_w;

%TE
rho = 1i.*k_0.*sigma_hat;
%TM
tau2 = epsilon_u./epsilon_w;
eta = sigma_hat./(1i.*k_0.*epsilon_w);

% g_bar = 0.025;
g_bar = 100;

Delta_TE = zeros(M,P+1);
Delta_TM = zeros(M,P+1);

for pp = 0:P
    Delta_TE(:,pp+1) = abs((g_bar.*k_w).*diff_besselj(pp,1,g_bar.*k_w)./besselj(pp,g_bar.*k_w)...
        -rho - (g_bar.*k_u).*diff_besselh(pp,1,g_bar.*k_u)./besselh(pp,g_bar.*k_u));
    Delta_TM(:,pp+1) = abs(tau2.*(g_bar.*k_w).*diff_besselj(pp,1,g_bar.*k_w)./besselj(pp,g_bar.*k_w)...
        + (g_bar.*k_u).*diff_besselh(pp,1,g_bar.*k_u)./besselh(pp,g_bar.*k_u)...
        .*(-1+eta.*(g_bar.*k_w).*diff_besselj(pp,1,g_bar.*k_w)./besselj(pp,g_bar.*k_w)));
   
end

% figure(1)    
% semilogy(lambda,Delta_TE(:,1),'b-o',lambda,Delta_TE(:,2),'g-*',lambda,Delta_TE(:,3),'m-x',...
%     lambda,Delta_TE(:,4),'r-s',lambda,Delta_TE(:,5),'y-d');
% xlabel('$\lambda$','interpreter','latex')
% ylabel('$\Delta_p$','interpreter','latex');
% title('$\Delta_p$(TE)','interpreter','latex');
% legend('$p=0$','$p=1$','$p=2$','$p=3$','$p=4$');
    

figure(2)    
semilogy(lambda,Delta_TM(:,1),'b-o',lambda,Delta_TM(:,2),'g-*',lambda,Delta_TM(:,3),'m-x',...
    lambda,Delta_TM(:,4),'r-s',lambda,Delta_TM(:,5),'y-d');
xlabel('$\lambda$','interpreter','latex')
ylabel('$\Delta_p$','interpreter','latex');
title('$\Delta_p$(TM)','interpreter','latex');
legend('$p=0$','$p=1$','$p=2$','$p=3$','$p=4$');

% figure(3)
% semilogy(lambda,Delta_TE(:,2),'b-o')
