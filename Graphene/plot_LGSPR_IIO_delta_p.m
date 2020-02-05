% PLOT_LGSPR_DELTA_P.m 
%
% Script to plot the determinant of LGSPR
%
% XT 11/19

close all

SavePlots = 0;
L = 2*pi;
N_theta = 64;
N_lambda = 101;
N_gbar = 100;
pp = 1;

lambda_low = 25; lambda_high = 50;
lambda = linspace(lambda_low,lambda_high,N_lambda)';
g_bar_low = 0.01; g_bar_high = 5;
g_bar = linspace(g_bar_low,g_bar_high,N_gbar)';

k_0 = L./lambda;
n_u = zeros(N_lambda,1); epsilon_u = zeros(N_lambda,1);
n_w = zeros(N_lambda,1); epsilon_w = zeros(N_lambda,1);
Delta_TM = zeros(N_lambda,N_gbar);
sigma_hat = zeros(N_lambda,1);
% IN = 'SILVER';
IN = 'VACUUM';
OUT = 'VACUUM';


omega = k_0*(3e8)*(1e6);
for m = 1:N_lambda
    [n_u(m),epsilon_u(m)] = ri_perm(lambda(m),OUT);
    [n_w(m),epsilon_w(m)] = ri_perm(lambda(m),IN);
%     [sigma_hat(m),~,~,~,~]...
%         = sigma_hat_graphene_low(omega(m),300,0.4,1,0.45,2.6*(1e-3));
%         %ALMA
    [~,~,~,~,sigma_hat(m)]...
        = sigma_hat_graphene_low(omega(m),300,0.4,1,0.45,2.6*(1e-3)); %BFPV
end
k_u = n_u.*k_0; 
k_w = n_w.*k_0;

%TM
tau2 = epsilon_u./epsilon_w;
eta = sigma_hat./(1i.*k_0.*epsilon_w);

for n = 1:N_gbar
    r = g_bar(n);
%     Delta_TM(:,n) = abs(tau2.*(r.*k_w).*diff_besselj(pp,1,r.*k_w)./besselj(pp,r.*k_w)...
%         + (r.*k_u).*diff_besselh(pp,1,r.*k_u)./besselh(pp,r.*k_u)...
%         .*(-1+eta.*(r.*k_w).*diff_besselj(pp,1,r.*k_w)./besselj(pp,r.*k_w)));
    Delta_TM(:,n) = abs(real(tau2.*(r.*k_w).*diff_besselj(pp,1,r.*k_w)./besselj(pp,r.*k_w)...
        + (r.*k_u).*diff_besselh(pp,1,r.*k_u)./besselh(pp,r.*k_u)...
        .*(-1+eta.*(r.*k_w).*diff_besselj(pp,1,r.*k_w)./besselj(pp,r.*k_w))));
end


figure(1)  
contourf(lambda,g_bar,log(Delta_TM).',40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\bar{g}$','interpreter','latex');
% set(gca,'ZScale', 'log')
title('LGSPR $\Delta_p$ versus $\lambda$ and $\bar{g}$','interpreter','latex');

% x1=0.5*(lambda_low+lambda_high); y1=0.4*(g_bar_low+g_bar_high);
% x2=0.5*(lambda_low+lambda_high); y2=0.6*(g_bar_low+g_bar_high);
% txt1 = ['lambda= ' num2str(lambda_low) ' to ' num2str(lambda_high)];
% txt2 = ['radius= ' num2str(g_bar_low) ' to ' num2str(g_bar_high) ];
% text(x1,y1,txt1)
% text(x2,y2,txt2)

if(SavePlots==1)
    filename = sprintf('fig_LGSPR_Delta_lambda%.0f_radius%.0f_%s%s',lambda_high,g_bar_high,IN,OUT);
    saveas(1,filename,'epsc');
end

