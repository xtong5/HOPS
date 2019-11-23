% test_dno_tfe_TE_twolayer.m
%
% Script to test TFE solvers for Graphene TE mode using DNO
%
% XT 8/19

clear all;
close all;

warning off;
SavePlots = 0;
SaveData = 0;

L = 2*pi;
N_theta = 64;

lambda = 0.45;
k_0 = L./lambda;
n_u = 1; %vacuum
n_w = 1; %vacuum
k_u = n_u.*k_0; 
k_w = n_w.*k_0;
epsilon_u = n_u.^2;
epsilon_w = n_w.^2;

[sigma_hat_alma,sigma_hat_zsdlz,sigma_hat_drude,...
    sigma_hat_kubo,sigma_hat_drude_bfpv]...
    = sigma_hat_graphene_low(k_0,0,0,0,0.45,2.6*(1e-3));
sigma_hat = sigma_hat_drude_bfpv;
rho = 1i.*k_0.*sigma_hat; %TE

a = 1.2; b = 1.6; c = 0.6;
N = 16;
% N = 0;
N_r = 16;
Eps = 0.02;
% Eps = 0;


fprintf('test_dno_tfe_TE_twolayer\n');
fprintf('---------------------------------------------\n');
fprintf('k_u = %g  k_w = %g sigma_hat = %g + %gi\n',...
    k_u,k_w,real(sigma_hat),imag(sigma_hat));
fprintf('Eps = %g  a = %g  b = %g  c = %g\n',Eps,a,b,c);
fprintf('N_theta = %d N = %d  N_r = %d\n',N_theta,N,N_r);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = ifft(1i*p.*fft(f));

A=a+Eps.*f;

Ar_u = 2; pp = 2; % take a special wavenumber
Ar_w = 1; %r = 2; % take a special wavenumber ???same

xi_u = Ar_u*besselh(pp,k_u.*A).*exp(1i*pp.*theta);
nu_u = Ar_u*((-k_u.*A.*(diff_besselh(pp,1,k_u.*A))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k_u.*A)./A ).*exp(1i*pp.*theta)); 
xi_w = Ar_w*besselj(pp,k_w.*A).*exp(1i*pp.*theta);
nu_w = Ar_w*((k_w.*A.*(diff_besselj(pp,1,k_w.*A))-...
    1i*pp*Eps.*f_theta.*besselj(pp,k_w.*A)./A ).*exp(1i*pp.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar_u*besselh(pp,k_u*a).*exp(1i*pp.*theta);
xi_u_n(:,1+1) = Ar_u*k_u^1*diff_besselh(pp,1,k_u*a).*f_n.*exp(1i*pp.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_besselh(pp,1,k_u*a).*exp(1i*pp.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_besselh(pp,1+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_besselh(pp,1,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar_u*(f_theta/a).*(1i*pp).*besselh(pp,k_u*a).*f_nmo.*exp(1i*pp.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(pp,k_w*a).*exp(1i*pp.*theta);
xi_w_n(:,1+1) = Ar_w*k_w^1*diff_besselj(pp,1,k_w*a).*f_n.*exp(1i*pp.*theta);
nu_w_n(:,0+1) = Ar_w*k_w*a*diff_besselj(pp,1,k_w*a).*exp(1i*pp.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_besselj(pp,1+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_besselj(pp,1,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar_w*(f_theta/a).*(1i*pp).*besselj(pp,k_w*a).*f_nmo.*exp(1i*pp.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_besselh(pp,n,k_u*a).*f_n.*exp(1i*pp.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_besselh(pp,n+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_besselh(pp,n,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_besselh(pp,n-1,k_u*a).*f_nmt.*exp(1i*pp.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*pp).*diff_besselh(pp,n-1,k_u*a)...
      .*f_nmo.*exp(1i*pp.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_besselj(pp,n,k_w*a).*f_n.*exp(1i*pp.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_besselj(pp,n+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_besselj(pp,n,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_besselj(pp,n-1,k_w*a).*f_nmt.*exp(1i*pp.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*pp).*diff_besselj(pp,n-1,k_w*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end

% BCs
zeta_n = xi_u_n - xi_w_n; % nu_u points downwards!
psi_n = -nu_u_n - nu_w_n + rho*current_n(f_theta,xi_w_n,N);


% Two-layer scattering by DNO


fprintf('\n\nTwo-layer scattering by DNO\n\n');

tic;
[U_n,W_n] = dno_tfe_TE_twolayer(zeta_n,psi_n,f,f_theta,rho,...
    p,k_u,k_w,a,b,c,N_theta,N,N_r);
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(U_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r);
Gn_tfe_u = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N,N_r);
[Wn,Dr_Wn,Dp_Wn] = field_tfe_helmholtz_polar_interior(W_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r);
Gn_tfe_w = dno_tfe_helmholtz_polar_interior(Dr_Wn,Dp_Wn,f,f_theta,k_w,a,c,p,N_theta,N,N_r);
t_tfe = toc;

% if SaveData==1
% % filename = sprintf('DNO_fe_Eps_%g.mat',Eps);
% filename = sprintf('DNO_fe_Eps_%g_sing12.mat',Eps);
% save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
%     'U_n','W_n','Gn_fe_u','Gn_fe_w','xi_u','xi_w','nu_u','nu_w');
% end

fprintf('Press key to compute exterior layer errors...\n');
pause;

fprintf('  t_tfe = %g\n',t_tfe);
fprintf('\nEXTERIOR LAYER\n\n');
[relerrU,nplotU] = compute_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerrDNOU,nplotDNOU] = compute_errors_2d_polar(nu_u,Gn_tfe_u,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotU,relerrU);
% make_plots_polar(SavePlots,nplotDNOU,relerrDNOU);
% fprintf('\n');
% 
fprintf('Press key to compute interior layer errors...\n');
pause;

fprintf('\nINTERIOR LAYER\n\n');
[relerrW,nplotW] = compute_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerrDNOW,nplotDNOW] = compute_errors_2d_polar(nu_w,Gn_tfe_w,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotW,relerrW);
% make_plots_polar(SavePlots,nplotDNOW,relerrDNOW);
% fprintf('\n');

