% PLOT_DELTA_P.m 
%
% Script to plot the determinant of 
%
% XT 9/19

close all

SavePlots = 0;
L = 2*pi;
N_theta = 64;
N_lambda = 1001;
N_gbar = 1001;
pp = 1;

lambda_low = 0.3; lambda_high = 0.8;
lambda = linspace(lambda_low,lambda_high,N_lambda)';
g_bar_low = 0.01; g_bar_high = 100;
g_bar = linspace(g_bar_low,g_bar_high,N_gbar)';

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

Delta_TM = zeros(N_lambda,N_gbar);

for n = 1:N_gbar
    r = g_bar(n);
    Delta_TM(:,n) = abs(-tau2.*(r.*k_w).*besselh(pp,r.*k_u).*diff_besselj(pp,1,r.*k_w)./besselj(pp,r.*k_w)...
        - (r.*k_u).*besselj(pp,r.*k_w).*diff_besselh(pp,1,r.*k_u)./besselh(pp,r.*k_u));
end

   

figure(1)  
contourf(lambda,g_bar,log(Delta_TM),40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\bar{g}$','interpreter','latex');
title('LSPR $\Delta_p$ versus $\lambda$ and $\bar{g}$','interpreter','latex');

% x1=0.5*(lambda_low+lambda_high); y1=0.4*(g_bar_low+g_bar_high);
% x2=0.5*(lambda_low+lambda_high); y2=0.6*(g_bar_low+g_bar_high);
% txt1 = ['lambda= ' num2str(lambda_low) ' to ' num2str(lambda_high)];
% txt2 = ['radius= ' num2str(g_bar_low) ' to ' num2str(g_bar_high) ];
% text(x1,y1,txt1)
% text(x2,y2,txt2)

if(SavePlots==1)
    filename = sprintf('fig_LSPR_Delta_lambda_radius%.0f_%s%s',g_bar_high,IN,OUT);
    saveas(1,filename,'epsc');
end

