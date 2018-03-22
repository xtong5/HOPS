clear all;
close all;
warning off;

%% Default
k=2;
Eps = 0.02;
% Eps = 0;
N_theta = 64;
N_r = 16;
a = 2;
c = 1;
N = 16; 
L = 2*pi;
p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';
fprintf('test_helmholtz_polar_TFE\n');
fprintf('-------------\n');
fprintf('k = %g a = %g c = %g\n',k,a,c);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');

%% pre allocation
u_n = zeros(N_theta*(N_r+1),N+1);
u_n_p = zeros(N_theta*(N_r+1),N+1);
Dr_u_n = zeros((N_r+1)*N_theta,N+1);
Dp_u_n = zeros((N_r+1)*N_theta,N+1);
J_n_p = zeros(N_theta,N+1);
Un = zeros(N_theta,N+1);
Dr_Un = zeros(N_theta,N+1);
Dp_Un = zeros(N_theta,N+1);
Gn_tfe = zeros(N_theta,N+1);

%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

%% Exact Solution
Ar = 1; pp = 2; % compute a special wavenumber
xi = Ar*besselj(pp,k.*(a+Eps.*f)).*exp(1i*pp.*theta);
u_exact = Ar * exp(1i*pp*theta) .* besselj(pp,k.*(a+Eps.*f));
f_n = ones(N_theta,1); 
xi_n(:,0+1) = Ar*besselj(pp,k*a).*exp(1i*pp.*theta);
for n=1:N
  f_n = f.*f_n/n;
  xi_n(:,n+1) = Ar*k^n*diff_besselj(pp,n,k*a).*f_n.*exp(1i*pp.*theta);
end
nu =Ar*(k.*(a+Eps.*f).*(diff_besselj(pp,1,k.*(a+Eps.*f)))-1i*pp*Eps.*...
    f_theta.*besselj(pp,k.*(a+Eps.*f))./(a+Eps.*f) ).*exp(1i*pp.*theta); %DNO

%% Main
tic;
dnp = field_fe_helmholtz_polar_interior(xi_n,f,k,a,p,N_theta,N);
% dnp1 = field_fe_helmholtz_polar_interior1(xi,f,k,p,N_theta,N,a);
Gn_fe = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k,a,p,N_theta,N);
t_fe = toc;

Un = field_tfe_helmholtz_polar_interior(xi_n,f,f_theta,k,a,c,p,N_theta,N,N_r);

Gn_tfe = dno_tfe_helmholtz_polar_interior(Un,f,f_theta,k,a,c,p,N_theta,N,N_r);
[relerr,nplot] = compute_errors_2d_polar(nu,Gn_tfe,Eps,N,N_theta); 
