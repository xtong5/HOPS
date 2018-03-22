% test_IIO_TFE_helmholtz_polar_interior.m
%
% Script to test Helmholtz IIO solvers in polar (interior)
%
% XT 3/18

clear all
%% Default
SavePlots = 0;
% RunNumber = 1;
Mode = 2; 
L = 2*pi;
k=2.1;
if(Mode==1)
  sigma = 1;
else
  sigma = (0.4*k/L)^2;
end
% sigma = Inf;
eta = 1.5;
% eta = -1i;

Eps = 0.02;
% Eps = 0;
N_theta = 64;
N = 16;
N_r = 16;
a = 1.6;
c = 0.9;


p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_IIO_TFE_helmholtz_polar_interior\n');
fprintf('-------------\n');
fprintf('k = %g a = %g c = %g  sigma = %g  eta = %g\n',k,a,c,sigma,eta);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');


%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
AA=a+Eps.*f;
u_exact = Ar * besselj(pp,k.*AA).*exp(1i*pp.*theta);
xi_w = Ar*besselj(pp,k.*AA).*exp(1i*pp.*theta);
nu_w = Ar*((k.*AA.*(diff_besselj(pp,1,k.*AA))-...
    1i*pp*Eps.*f_theta.*besselj(pp,k.*AA)./AA ).*exp(1i*pp.*theta)); 
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_w_n(:,0+1) = Ar*besselj(pp,k*a).*exp(1i*pp.*theta);
xi_w_n(:,1+1) = Ar*k^1*diff_besselj(pp,1,k*a).*f_n.*exp(1i*pp.*theta);
nu_w_n(:,0+1) = Ar*k*a*diff_besselj(pp,1,k*a).*exp(1i*pp.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar*a*k^(1+1).*diff_besselj(pp,1+1,k*a).*f_n.*exp(1i*pp.*theta)...
      +Ar*(2*f).*k^1.*diff_besselj(pp,1,k*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar*(f_theta/a).*(1i*pp).*besselj(pp,k*a).*f_nmo.*exp(1i*pp.*theta);
for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_w_n(:,n+1) = Ar*k^n*diff_besselj(pp,n,k*a).*f_n.*exp(1i*pp.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar*a*k^(n+1).*diff_besselj(pp,n+1,k*a).*f_n.*exp(1i*pp.*theta)...
      +Ar*(2*f).*k^n.*diff_besselj(pp,n,k*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar*(f.^2/a)*k^(n-1).*diff_besselj(pp,n-1,k*a).*f_nmt.*exp(1i*pp.*theta)...
      -Ar*(f_theta/a)*k^(n-1).*(1i*pp).*diff_besselj(pp,n-1,k*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end

I_w = 1/sigma*nu_w+1i*eta*xi_w;
I_w_n = 1/sigma*nu_w_n+1i*eta*xi_w_n;
S_w = 1/sigma*nu_w-1i*eta*xi_w;
Sn_w = 1/sigma*nu_w_n-1i*eta*xi_w_n;

% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k,a,c,p,N_theta,N,N_r,sigma,eta);
% [relerr,nplot] = compute_errors_2d_polar(u_exact,Un,Eps,N,N_theta);
% 
% Sn_tfe = IIO_tfe_helmholtz_polar_interior(Un,Dr_Un,Dp_Un,f,f_theta,a,c,N,sigma,eta);
% [relerr,nplot] = compute_errors_2d_polar(S_u,Sn_tfe,Eps,N,N_theta);        
% 

Sn_tfe = IIO_tfe_helmholtz_polar_interior1(I_w_n,f,f_theta,k,a,c,p,N_theta,N,N_r,sigma,eta);
[relerr,nplot] = compute_errors_2d_polar(S_w,Sn_tfe,Eps,N,N_theta);        
