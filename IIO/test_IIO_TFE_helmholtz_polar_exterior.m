% test_IIO_TFE_helmholtz_polar_exterior.m
%
% Script to test Helmholtz IIO solvers in polar (exterior)
%
% XT 3/18


%% Default
SavePlots = 0;
% RunNumber = 1;
Mode = 2; 
L = 2*pi;
lambda = 0.4;
% k=2.1;
k=13.962634015954636;
if(Mode==1)
  sigma = 1;
else
  sigma = (lambda*k/L)^2;
end
eta = 3.4;

Eps = 0.02;
% Eps = 0;
N_theta = 64;
N = 16;
N_r = 16;
a = 1.2;
b = 1.6;


p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_IIO_TFE_helmholtz_polar_exterior\n');
fprintf('-------------\n');
fprintf('k = %g a = %g b = %g  sigma = %g  eta = %g\n',k,a,b,sigma,eta);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');


%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
AA=a+Eps.*f;
% xi = exp(1i*pp.*theta).*(-1/sigma*(k.*AA*Ar.*diff_besselh(pp,1,k.*AA)...
%     -Eps.*f_theta./AA*(1i*pp)*Ar.*besselh(pp,k.*AA))...
%     +1i*eta*Ar.*besselh(pp,k.*AA));
u_exact = Ar * exp(1i*pp*theta) .* besselh(pp,k.*(a+Eps.*f));
xi_u = Ar*besselh(pp,k.*AA).*exp(1i*pp.*theta);
nu_u = Ar*((-k.*AA.*(diff_besselh(pp,1,k.*AA))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k.*AA)./AA ).*exp(1i*pp.*theta)); 
I_u = 1/sigma*nu_u+1i*eta*xi_u;
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar*besselh(pp,k*a).*exp(1i*pp.*theta);
xi_u_n(:,1+1) = Ar*k^1*diff_besselh(pp,1,k*a).*f_n.*exp(1i*pp.*theta);
nu_u_n(:,0+1) = -Ar*k*a*diff_besselh(pp,1,k*a).*exp(1i*pp.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar*a*k^(1+1).*diff_besselh(pp,1+1,k*a).*f_n.*exp(1i*pp.*theta)...
      -Ar*(2*f).*k^1.*diff_besselh(pp,1,k*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar*(f_theta/a).*(1i*pp).*besselh(pp,k*a).*f_nmo.*exp(1i*pp.*theta);
for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar*k^n*diff_besselh(pp,n,k*a).*f_n.*exp(1i*pp.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar*a*k^(n+1).*diff_besselh(pp,n+1,k*a).*f_n.*exp(1i*pp.*theta)...
      -Ar*(2*f).*k^n.*diff_besselh(pp,n,k*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar*(f.^2/a)*k^(n-1).*diff_besselh(pp,n-1,k*a).*f_nmt.*exp(1i*pp.*theta)...
      +Ar*(f_theta/a)*k^(n-1).*(1i*pp).*diff_besselh(pp,n-1,k*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end

I_u_n = 1/sigma*nu_u_n+1i*eta*xi_u_n;
Q_u = 1/sigma*nu_u-1i*eta*xi_u;
Qn_u = 1/sigma*nu_u_n-1i*eta*xi_u_n;

[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k,a,b,p,N_theta,N,N_r,sigma,eta);
[relerr,nplot] = compute_errors_2d_polar(u_exact,Un,Eps,N,N_theta);
 
Qn_tfe = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma,eta);
[relerr,nplot] = compute_errors_2d_polar(Q_u,Qn_tfe,Eps,N,N_theta);        



