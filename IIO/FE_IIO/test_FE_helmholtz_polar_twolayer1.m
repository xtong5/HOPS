% test_helmholtz_polar_twolayer1.m%
%
% Script to find a bad example (two layers)
%
% XT 11/17

clear all;
close all;
SavePlots = 0;

RunNumber = 3;
Mode = 2; %check 
% Mode = 1;

L = 2*pi;
lambda = 0.45;
n_u = 1;
n_w = 2.5;
k_0 = L/lambda;
k_u = n_u*k_0; 
k_w = n_w*k_0;

if(RunNumber==1)
  % Small Deformation
  Eps = 0.002;
  N_theta = 64;
  a = 0.025;
  b = 10*a;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  a = 2.0;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.02;
  N_theta = 64;
  a = 0.025;
  N = 24;
end

fprintf('test_helmholtz_twolayer_polar\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;
A=a+Eps.*f;

Ar_u = 2; r = 2; % take a special wavenumber
Ar_w = 1; %r = 2; % take a special wavenumber ???same

xi_u = Ar_u*besselh(r,k_u.*A).*exp(1i*r.*theta);
nu_u = Ar_u*((-k_u.*A.*(diff_bessel(2,r,1,k_u.*A))+...
    1i*r*Eps.*f_theta.*besselh(r,k_u.*A)./A ).*exp(1i*r.*theta)); 
xi_w = Ar_w*besselj(r,k_w.*A).*exp(1i*r.*theta);
nu_w = Ar_w*((k_w.*A.*(diff_bessel(1,r,1,k_w.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k_w.*A)./A ).*exp(1i*r.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar_u*besselh(r,k_u*a).*exp(1i*r.*theta);
xi_u_n(:,1+1) = Ar_u*k_u^1*diff_bessel(2,r,1,k_u*a).*f_n.*exp(1i*r.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_bessel(2,r,1,k_u*a).*exp(1i*r.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_bessel(2,r,1+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_bessel(2,r,1,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a).*(1i*r).*besselh(r,k_u*a).*f_nmo.*exp(1i*r.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(r,k_w*a).*exp(1i*r.*theta);
xi_w_n(:,1+1) = Ar_w*k_w^1*diff_bessel(1,r,1,k_w*a).*f_n.*exp(1i*r.*theta);
nu_w_n(:,0+1) = Ar_w*k_w*a*diff_bessel(1,r,1,k_w*a).*exp(1i*r.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_bessel(1,r,1+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_bessel(1,r,1,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a).*(1i*r).*besselj(r,k_w*a).*f_nmo.*exp(1i*r.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_bessel(2,r,n,k_u*a).*f_n.*exp(1i*r.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_bessel(2,r,n+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_bessel(2,r,n,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_bessel(2,r,n-1,k_u*a).*f_nmt.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*r).*diff_bessel(2,r,n-1,k_u*a)...
      .*f_nmo.*exp(1i*r.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_bessel(1,r,n,k_w*a).*f_n.*exp(1i*r.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_bessel(1,r,n+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_bessel(1,r,n,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_bessel(1,r,n-1,k_w*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*r).*diff_bessel(1,r,n-1,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);
end

% testing
%xi_u_n = zeros(N_theta,N+1);
%nu_u_n = zeros(N_theta,N+1);
%xi_w_n = zeros(N_theta,N+1);
%nu_w_n = zeros(N_theta,N+1);
%xi_u_n(:,0+1) = xi_u;
%nu_u_n(:,0+1) = nu_u;
%xi_w_n(:,0+1) = xi_w;
%nu_w_n(:,0+1) = nu_w;
% end testing

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;
%psi_n = nu_u_n - tau2*nu_w_n;


% Two-layer scattering by DNO


fprintf('\n\nTwo-layer scattering by DNO\n\n');

tic;
U_n = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N);
apn_fe = field_fe_helmholtz_polar_exterior(U_n,f,k_u,a,p,N_theta,N);
Gn_fe_u = dno_fe_helmholtz_polar_exterior(apn_fe,f,f_theta,k_u,a,p,N_theta,N);
W_n = U_n - zeta_n;
dpn_fe = field_fe_helmholtz_polar_interior(W_n,f,k_w,a,p,N_theta,N);
Gn_fe_w = dno_fe_helmholtz_polar_interior(dpn_fe,f,f_theta,k_w,a,p,N_theta,N);
t_fe = toc;

fprintf('Press key to compute exterior layer errors...\n');

fprintf('  t_fe = %g\n',t_fe);
fprintf('\nEXTERIOR LAYER\n\n');
[relerrDNOU,nplotDNOU] = compute_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
make_plots_polar(SavePlots,nplotDNOU,relerrDNOU);
fprintf('\n');

% filename = sprintf('FE_Eps_%g_N%g.mat',Eps,N);
% save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a','Gn_fe_u');
% 
% 

