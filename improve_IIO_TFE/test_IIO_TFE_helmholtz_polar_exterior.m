% test_IIO_TFE_helmholtz_polar_exterior.m
%
% Script to test Helmholtz IIO solvers in polar (exterior)
%
% XT 3/18 5/19


%% Default
SavePlots = 0;
Mode = 1; 
L = 2*pi;
lambda = 0.4;
n_u = 1.2;
k_0 = L/lambda;
k_u = n_u*k_0; 

if(Mode==1)
  sigma_u = 1;
else
  sigma_u = 1./(lambda*k_u/L)^2;
end
eta = 3.4;

N_theta = 64;
N = 16;
N_r = 16;
% a = 0.5;b = 0.8;
% Eps = 0.02;
a = 0.025;b = 2*a;%10a
Eps = 0.01*a;
% Eps = 0;


p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_IIO_TFE_helmholtz_polar_exterior\n');
fprintf('-------------\n');
fprintf('k = %g a = %g b = %g  sigma = %g  eta = %g\n',k_u,a,b,sigma_u,eta);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
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

Y_p = 1i*eta*ones(N_theta,1);
Z_p = -1i*eta*ones(N_theta,1);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_u_n = sigma_u*nu_u_n+ifft(Y_p.*fft(xi_u_n));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
Qn_u = sigma_u*nu_u_n+ifft(Z_p.*fft(xi_u_n));

tic
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,Y_p);
[relerr,nplot] = compute_errors_2d_polar(u_exact,Un,Eps,N,N_theta);
 
Qn_tfe = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,Z_p);
[relerr,nplot] = compute_errors_2d_polar(Q_u,Qn_tfe,Eps,N,N_theta);        
t = toc;

fprintf('  t_tfe = %g\n',t);


