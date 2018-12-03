% test_helmholtz_polar_exterior.m
%
% Script to test Helmholtz IIO FE solvers in polar (exterior)
%
% XT 10/18

% clear all;
% close all;

%% Default
SavePlots = 0;
% RunNumber = 1;
warning off 

Mode = 2; 
L = 2*pi;
k_u=2.1;
% k_u = 1;
k_w=1.1;
if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = (0.4*k_u/L)^2;
  sigma_w = (0.4*k_w/L)^2;
end

Eps = 0.02;
% Eps = 0;
N_theta = 64;
N = 16;
% N = 0;
a = 2.1;
% a = 1;
% sigma_u = inf;

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_IIO_TFE_helmholtz_polar_exterior\n');
fprintf('-------------\n');
fprintf('k_u = %g a = %g  sigma = %g\n',k_u,a,sigma_u);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');


%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
AA=a+Eps.*f;
u_exact = Ar * exp(1i*pp*theta) .* besselh(pp,k_u.*(a+Eps.*f));
xi_u = Ar*besselh(pp,k_u.*AA).*exp(1i*pp.*theta);
nu_u = Ar*((-k_u.*AA.*(diff_besselh(pp,1,k_u.*AA))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k_u.*AA)./AA ).*exp(1i*pp.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar*besselh(pp,k_u*a).*exp(1i*pp.*theta);
xi_u_n(:,1+1) = Ar*k_u^1*diff_besselh(pp,1,k_u*a).*f_n.*exp(1i*pp.*theta);
nu_u_n(:,0+1) = -Ar*k_u*a*diff_besselh(pp,1,k_u*a).*exp(1i*pp.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar*a*k_u^(1+1).*diff_besselh(pp,1+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar*(2*f).*k_u^1.*diff_besselh(pp,1,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar*(f_theta/a).*(1i*pp).*besselh(pp,k_u*a).*f_nmo.*exp(1i*pp.*theta);
for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar*k_u^n*diff_besselh(pp,n,k_u*a).*f_n.*exp(1i*pp.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar*a*k_u^(n+1).*diff_besselh(pp,n+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar*(2*f).*k_u^n.*diff_besselh(pp,n,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar*(f.^2/a)*k_u^(n-1).*diff_besselh(pp,n-1,k_u*a).*f_nmt.*exp(1i*pp.*theta)...
      +Ar*(f_theta/a)*k_u^(n-1).*(1i*pp).*diff_besselh(pp,n-1,k_u*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end

Z_p = k_u * diff_besselh(p,1,k_u*a)./(sigma_u*besselh(p,k_u*a));
Y_p = k_w * diff_besselj(p,1,k_w*a)./(sigma_w*besselj(p,k_w*a));

% Y_p = 1i*3.4.*ones(N_theta,1);
% Y_p = 1.2*ones(N_theta,1);
% Z_p = 1.1*ones(N_theta,1);


I_u = 1/sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_u_n = 1/sigma_u*nu_u_n+ifft(Y_p.*fft(xi_u_n));
Q_u = 1/sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
% Qn_u = 1/sigma_u*nu_u_n+ifft(Z_p.*fft(xi_u_n));


tic;
anp= field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Qn = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
t_fe = toc;
fprintf('  t_fe = %g\n',t_fe);
[relerr,nplot] = compute_errors_2d_polar(Q_u,Qn,Eps,N,N_theta);

% make_plots_polar(SavePlots,nplot,relerr);