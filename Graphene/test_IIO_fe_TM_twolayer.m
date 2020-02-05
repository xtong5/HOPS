% test_dno_fe_TM_twolayer.m
%
% Script to test IIO solvers in polar (two layers) for Graphene TM mode 
%
% XT 9/19

clear all;
close all;

warning off;
SavePlots = 0;
SaveData = 1;

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
tau2 = epsilon_u./epsilon_w;


%TM
[~,~,~,~,sigma_hat]...
    = sigma_hat_graphene_low(k_0,0,0,0,0.45,2.6*(1e-3));
eta = sigma_hat./(1i.*k_0.*epsilon_w);
sigma_u = 1./epsilon_u; sigma_w = 1./epsilon_w;
% sigma_u = 1; sigma_w = 1;

Y_p = 1i*3.4.*ones(N_theta,1); Z_p = -1i*3.4.*ones(N_theta,1);
a = 0.025;
Eps = 0.01*a;
N = 16;
% Eps = 0;


fprintf('test_IIO_TM_twolayer\n');
fprintf('---------------------------------------------\n');
fprintf('k_u = %g  k_w = %g sigma_hat = %g\n',k_u,k_w,sigma_hat);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d\n',N_theta,N);
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
Nzeta_n = current_n(f_theta,xi_u_n,N)-current_n(f_theta,xi_w_n,N) + eta.*nu_w_n;
psi_n = -nu_u_n - tau2.*nu_w_n; % nu_u points downwards!

% exact 
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));


% Two-layer scattering by IIO

fprintf('\n\nTwo-layer scattering by IIO\n\n');

tic;
[I_u_n,I_w_n] = IIO_fe_TM_twolayer(Nzeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,eta,a,N_theta,N,Y_p,Z_p);
anp = field_fe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

if SaveData==1
filename = sprintf('IIO_fe_TM_Eps_%g.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','a',...
    'I_u_n','I_w_n','Q_u_n','S_w_n','I_u','I_w','Q_u','S_w');
end

fprintf('Press key to compute exterior layer errors...\n');
pause;

fprintf('  t_fe = %g\n',t_fe);
fprintf('\nEXTERIOR LAYER\n\n');
[relerrIU,nplotIU] = compute_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerrIIOU,nplotIIOU] = compute_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotU,relerrU);
% make_plots_polar(SavePlots,nplotIIOU,relerrIIOU);
fprintf('\n');

fprintf('Press key to compute interior layer errors...\n');
pause;

fprintf('\nINTERIOR LAYER\n\n');
[relerrIW,nplotIW] = compute_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerrIIOW,nplotIIOW] = compute_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotW,relerrW);
% make_plots_polar(SavePlots,nplotIIOW,relerrIIOW);
fprintf('\n');

