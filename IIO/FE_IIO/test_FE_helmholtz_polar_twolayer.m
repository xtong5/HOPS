% test_helmholtz_polar_twolayer.m
%
% Script to test Helmholtz IIO solvers in polar (two layers)
%
% XT 11/18

clear all;
close all;
warning off;
SavePlots = 0;

RunNumber = 1;
Mode = 1; %check 

L = 2*pi;
lambda = 0.45;
% n_u = 1;
% n_w = 2.5;
k_0 = L/lambda;
OUT = 'VACUUM';IN = 'SILVER';
[n_u,epsilon_u] = ri_perm(lambda,OUT);
[n_w,epsilon_w] = ri_perm(lambda,IN);
k_u = n_u*k_0; 
k_w = n_w*k_0;
% k_w = 5.13562230184068;
if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = 1/epsilon_u;
  sigma_w = 1/epsilon_w;
end

if(RunNumber==1)
  % Small Deformation
%   Eps = 0.02;
%   Eps = 0;
  N_theta = 64;
%   a = 1-1e-16; 
%   a = 1-1e-12; 
%   a = 0.5; 
  a = 0.025;
  Eps = 0.01*a;
  b = 10*a;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  a = 0.025;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  N_theta = 64;
  a = 0.025;
  N = 16;
end

fprintf('test_IIO_FE_twolayer\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;
AA=a+Eps.*f;

Ar_u = 2; Ar_w = 1; pp = 2; % take a special wavenumber
 
xi_u = Ar_u*besselh(pp,k_u.*AA).*exp(1i*pp.*theta);
nu_u = Ar_u*((-k_u.*AA.*(diff_besselh(pp,1,k_u.*AA))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k_u.*AA)./AA ).*exp(1i*pp.*theta)); 
xi_w = Ar_w*besselj(pp,k_w.*AA).*exp(1i*pp.*theta);
nu_w = Ar_w*((k_w.*AA.*(diff_besselj(pp,1,k_w.*AA))-...
    1i*pp*Eps.*f_theta.*besselj(pp,k_w.*AA)./AA ).*exp(1i*pp.*theta)); 
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

zeta_n = xi_u_n - xi_w_n; % nu_u points downwards!
psi_n = -nu_u_n.*sigma_u - nu_w_n.*sigma_w;

% Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
% Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

Z_p = -1i*3.4.*ones(N_theta,1);
Y_p = 1i*3.4.*ones(N_theta,1);


I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));
U_far = Ar_u*besselh(pp,k_u.*b).*exp(1i*pp.*theta); % true farfield 

% Two-layer scattering by IIO FE


fprintf('\n\nTwo-layer scattering by IIO\n\n');

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
Un_far = IIO_fe_farfield(anp,k_u,a,b,p,N_theta,N);
Un = IIO_fe_Un(anp,f,k_u,a,p,N_theta,N);
Wn = Un - zeta_n;
t_fe = toc;

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
% fprintf('\n');

fprintf('Press key to compute Dirichlet data errors...\n');
pause;

fprintf('\Interface\n\n');
[relerrU,nplotU] = compute_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerrW,nplotW] = compute_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);

fprintf('Press key to compute far-field errors...\n');
pause;

fprintf('\Far-field\n\n');
[relerrBU,nplotBU] = compute_errors_2d_polar(U_far,Un_far,Eps,N,N_theta);




% filename = sprintf('IIO_new_Eps_%g_Nr%g_sing16.mat',Eps,N_r);
% filename = sprintf('IIO_new_Eps_%g_Nr%g.mat',Eps,N_r);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','Z_p','Y_p','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w');




