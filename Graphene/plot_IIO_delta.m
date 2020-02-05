% PLOT_IIO_DELTA.m 
%
% Script to plot the determinant of M_TE and M_TM
%
% XT 11/19

close all
clear all

SavePlots = 1;
L = 2*pi;
N_theta = 64;
N_lambda = 2001;
P = 4;

lambda_low = 10;
lambda_high = 40;
g_bar = 0.025;

lambda = linspace(lambda_low,lambda_high,N_lambda)';
sigma_hat = zeros(N_lambda,1);

k_0 = L./lambda;
n_u = zeros(N_lambda,1); epsilon_u = zeros(N_lambda,1);
n_w = zeros(N_lambda,1); epsilon_w = zeros(N_lambda,1);
% IN = 'SILVER';
IN = 'VACUUM';
OUT = 'VACUUM';

omega = k_0*(3e8)*(1e6);
for m = 1:N_lambda
    [n_u(m),epsilon_u(m)] = ri_perm(lambda(m),OUT);
    [n_w(m),epsilon_w(m)] = ri_perm(lambda(m),IN);
%     [sigma_hat(m),~,~,~,~]...
%         = sigma_hat_graphene_low(omega(m),300,0.4,1,0.45,2.6*(1e-3)); %ALMA
    [~,~,~,~,sigma_hat(m)]...
        = sigma_hat_graphene_low(omega(m),300,0.4,1,0.45,2.6*(1e-3)); %BFPV
end
k_u = n_u.*k_0; 
k_w = n_w.*k_0;

sigma_u = 1./epsilon_u; sigma_w = 1./epsilon_w;
eta = 3.4;
Y_p = 1i*eta; Z_p = -1i*eta;

%TE
C1 = 1i.*k_0.*sigma_hat./(2*1i*eta);
%TM
C2 = 1i.*eta.*sigma_hat./(2*1i.*k_0);

Delta_TE = zeros(N_lambda,P+1);
Delta_TM = zeros(N_lambda,P+1);

for pp = 0:P
    A = -sigma_u.*g_bar.*k_u.*diff_besselh(pp,1,k_u.*g_bar)./besselh(pp,k_u.*g_bar);
    B = sigma_w.*g_bar.*k_w.*diff_besselj(pp,1,k_w.*g_bar)./besselj(pp,k_w.*g_bar);
    Q0_p = (A + Z_p)./(A + Y_p); S0_p = (B - Y_p)./(B - Z_p); 
    Delta_TE(:,pp+1) = abs(1-C1.*(1-S0_p)-Q0_p.*S0_p-C1.*Q0_p.*(1-S0_p));
    Delta_TM(:,pp+1) = abs(1-C2.*(1+S0_p)-Q0_p.*S0_p-C2.*Q0_p.*(1+S0_p));
   
end

figure(1)    
semilogy(lambda,Delta_TE(:,1),'b-o',lambda,Delta_TE(:,2),'g-*',lambda,Delta_TE(:,3),'m-x',...
    lambda,Delta_TE(:,4),'r-s',lambda,Delta_TE(:,5),'y-d');
xlabel('$\lambda$','interpreter','latex','FontSize',20);
ylabel('$\Delta_p$','interpreter','latex','FontSize',20);
title('$\Delta_p$(TE)','interpreter','latex','FontSize',20);
legend('$p=0$','$p=1$','$p=2$','$p=3$','$p=4$','interpreter','latex','FontSize',14);
    

figure(2)    
semilogy(lambda,Delta_TM(:,1),'b-o',lambda,Delta_TM(:,2),'g-*',lambda,Delta_TM(:,3),'m-x',...
    lambda,Delta_TM(:,4),'r-s',lambda,Delta_TM(:,5),'y-d');
xlabel('$\lambda$','interpreter','latex','FontSize',20);
ylabel('$\Delta_p$','interpreter','latex','FontSize',20);
title('$\Delta_p$(TM)','interpreter','latex','FontSize',20);
legend('$p=0$','$p=1$','$p=2$','$p=3$','$p=4$','interpreter','latex','FontSize',14);

% [min, index] = min(Delta_TM(:,2));
% fprintf('smallest value of delta_1 = %f at labmda = %f\n',min,lambda(index));

if(SavePlots==1)
%     filename = sprintf('fig_LGSPR_IIO_TE_Delta_radius%.0f_%s%s_ALMA',g_bar*1000,IN,OUT);
%     saveas(1,filename,'epsc');
%     filename = sprintf('fig_LGSPR_IIO_TM_Delta_radius%.0f_%s%s_ALMA',g_bar*1000,IN,OUT);
%     saveas(2,filename,'epsc');
    filename = sprintf('fig_LGSPR_IIO_TE_Delta_radius%.0f_%s%s_BFPV',g_bar*1000,IN,OUT);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_LGSPR_IIO_TM_Delta_radius%.0f_%s%s_BFPV',g_bar*1000,IN,OUT);
    saveas(2,filename,'epsc');
end

