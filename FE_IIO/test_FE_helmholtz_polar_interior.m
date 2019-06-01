% test_helmholtz_polar_interior.m
%
% Script to test Helmholtz IIO FE solvers in polar (interior)
%
% XT 11/18

% clear all;
% close all;

%% Default
SavePlots = 0;
% RunNumber = 1;
warning off 

Mode = 1; 
L = 2*pi;
lambda = 0.4;
n_u = 1;
n_w = 2.5;
k_0 = L/lambda;
k_u = n_u*k_0; 
% k_w = n_w*k_0;
k_w = 5.13562230184068;
if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = 1/(lambda*k_u/L)^2;
  sigma_w = 1/(lambda*k_w/L)^2;
end

Eps = 0.02;
% Eps = 0;
N_theta = 64;
N = 16;
% N = 0;
a = 1-1e-12;
% a = 1;
% sigma_u = inf;

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_helmholtz_polar_interior\n');
fprintf('-------------\n');
fprintf('k_w = %g a = %g  sigma = %g\n',k_u,a,sigma_w);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

%% data
f = exp(cos(theta));
f_theta = -sin(theta).*f;


%% Exact Solution
Ar = 1; r = 2; % compute a special wavenumber
A=a+Eps.*f;
xi_w = Ar*besselj(r,k_w.*A).*exp(1i*r.*theta);
nu_w =Ar*((k_w.*A.*(diff_besselj(r,1,k_w.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k_w.*A)./A ).*exp(1i*r.*theta)); %DNO
xi_w_n = zeros(N_theta,N); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
xi_w_n(:,0+1) = Ar*besselj(r,k_w*a).*exp(1i*r.*theta);
nu_w_n(:,0+1) = Ar*k_w*a*diff_besselj(r,1,k_w*a).*exp(1i*r.*theta);
f_n = f.*f_n;
xi_w_n(:,1+1) = Ar*k_w^1*diff_besselj(r,1,k_w*a).*f_n.*exp(1i*r.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar*a*k_w^(1+1).*diff_besselj(r,1+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar*(2*f).*k_w^1.*diff_besselj(r,1,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar*(f_theta/a).*(1i*r).*besselj(r,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);


for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_w_n(:,n+1) = Ar*k_w^n*diff_besselj(r,n,k_w*a).*f_n.*exp(1i*r.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar*a*k_w^(n+1).*diff_besselj(r,n+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar*(2*f).*k_w^n.*diff_besselj(r,n,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar*(f.^2/a)*k_w^(n-1).*diff_besselj(r,n-1,k_w*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar*(f_theta/a)*k_w^(n-1).*(1i*r).*diff_besselj(r,n-1,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);
end

Z_p = k_u * diff_besselh(p,1,k_u*a)./(sigma_u*besselh(p,k_u*a));
Y_p = k_w * diff_besselj(p,1,k_w*a)./(sigma_w*besselj(p,k_w*a));

% Y_p = 1i*3.4.*ones(N_theta,1);
% Z_p = 1i*3.4.*ones(N_theta,1);

I_w = 1/sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
I_w_n = 1/sigma_w*nu_w_n-ifft(Z_p.*fft(xi_w_n));
S_w = 1/sigma_w*nu_w-ifft(Y_p.*fft(xi_w));
% Sn_w = 1/sigma_w*nu_w_n-ifft(Y_p.*fft(xi_w_n));


tic;
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
Sn = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;
fprintf('  t_fe = %g\n',t_fe);
[relerr,nplot] = compute_errors_2d_polar(S_w,Sn,Eps,N,N_theta);

% make_plots_polar(SavePlots,nplot,relerr);