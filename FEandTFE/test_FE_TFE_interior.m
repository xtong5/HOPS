clear all;
close all;
warning off;

%% Default
SavePlots = 0;
RunNumber = 1;
L = 2*pi;
k=2.1;

if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  c = 0.6;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  c = 0.6;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  c = 0.6;
end

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_FE_TFE_interior\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k = %g a = %g c = %g\n',k,a,c);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');

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
Gn_fe = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k,a,p,N_theta,N);
t_fe = toc;
tic;
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_interior(xi_n,f,f_theta,k,a,c,p,N_theta,N,N_r);
Gn_tfe = dno_tfe_helmholtz_polar_interior(Dr_Un,Dp_Un,f,f_theta,k,a,c,p,N_theta,N,N_r);
t_tfe = toc;
% fprintf('  t_fe = %g  t_tfe = %g\n',t_fe,t_tfe);
[relerr,nplot] = compute_errors_2d_polar(nu,Gn_fe,Gn_tfe,Eps,N,N_theta);
make_plots_polar(SavePlots,nplot,relerr);
