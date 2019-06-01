% test_IIO_TFE_helmholtz_polar_interior.m
%
% Script to test Helmholtz IIO solvers in polar (interior)
%
% XT 3/18 5/19

clear all
%% Default
SavePlots = 0;
Mode = 2; 
L = 2*pi;
lambda = 0.4;
k_w=2.1;
if(Mode==1)
  sigma_w = 1;
else
  sigma_w = 1./(lambda*k_w/L)^2;
end
eta = 3.4;

N_theta = 64;
N = 16;
N_r = 16;
% a = 0.8;c = 0.4;
% Eps = 0.02;
a = 0.025;c = 0.1*a;
Eps = 0.01*a;
% Eps = 0;


p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_IIO_TFE_helmholtz_polar_interior\n');
fprintf('-------------\n');
fprintf('k = %g a = %g c = %g  sigma = %g  eta = %g\n',k_w,a,c,sigma_w,eta);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');


%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
AA=a+Eps.*f;
u_exact = Ar * besselj(pp,k_w.*AA).*exp(1i*pp.*theta);
xi_w = Ar*besselj(pp,k_w.*AA).*exp(1i*pp.*theta);
nu_w = Ar*((k_w.*AA.*(diff_besselj(pp,1,k_w.*AA))-...
    1i*pp*Eps.*f_theta.*besselj(pp,k_w.*AA)./AA ).*exp(1i*pp.*theta)); 
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_w_n(:,0+1) = Ar*besselj(pp,k_w*a).*exp(1i*pp.*theta);
xi_w_n(:,1+1) = Ar*k_w^1*diff_besselj(pp,1,k_w*a).*f_n.*exp(1i*pp.*theta);
nu_w_n(:,0+1) = Ar*k_w*a*diff_besselj(pp,1,k_w*a).*exp(1i*pp.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar*a*k_w^(1+1).*diff_besselj(pp,1+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar*(2*f).*k_w^1.*diff_besselj(pp,1,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar*(f_theta/a).*(1i*pp).*besselj(pp,k_w*a).*f_nmo.*exp(1i*pp.*theta);
for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_w_n(:,n+1) = Ar*k_w^n*diff_besselj(pp,n,k_w*a).*f_n.*exp(1i*pp.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar*a*k_w^(n+1).*diff_besselj(pp,n+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar*(2*f).*k_w^n.*diff_besselj(pp,n,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar*(f.^2/a)*k_w^(n-1).*diff_besselj(pp,n-1,k_w*a).*f_nmt.*exp(1i*pp.*theta)...
      -Ar*(f_theta/a)*k_w^(n-1).*(1i*pp).*diff_besselj(pp,n-1,k_w*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end

Y_p = 1i*eta*ones(N_theta,1);
Z_p = -1i*eta*ones(N_theta,1);

I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
I_w_n = sigma_w*nu_w_n-ifft(Z_p.*fft(xi_w_n));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));
Sn_w = sigma_w*nu_w_n-ifft(Y_p.*fft(xi_w_n));

tic
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,Z_p);
[relerr,nplot] = compute_errors_2d_polar(u_exact,Un,Eps,N,N_theta);

Sn_tfe = IIO_tfe_helmholtz_polar_interior(Un,Dr_Un,Dp_Un,f,f_theta,a,c,N,sigma_w,Y_p);
[relerr,nplot] = compute_errors_2d_polar(S_w,Sn_tfe,Eps,N,N_theta);        
t = toc;

fprintf('  t_tfe = %g\n',t);