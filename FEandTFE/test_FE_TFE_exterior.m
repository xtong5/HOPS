clear all

%% Default
SavePlots = 0;
RunNumber = 3;
L = 2*pi;
k=2.1;

if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  b = 1.6;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  b = 1.6;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  %Eps = 0.75;
  Eps = 1.0;
  N_theta = 64;
  %N = 16;
  %N_r = 16;
  N = 32;
  N_r = 32;
  a = 1;
  b = 2.0;
end


fprintf('test_FE_TFE_exterior\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k = %g a = %g b = %g\n',k,a,b);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

Ar = 2; r = 2; % compute a special wavenumber
A=a+Eps.*f;
xi = Ar*besselh(r,k.*A).*exp(1i*r.*theta);
nu =Ar*((-k.*A.*(diff_bessel(2,r,1,k.*A))+...
    1i*r*Eps.*f_theta.*besselh(r,k.*A)./A ).*exp(1i*r.*theta)); %DNO
xi_n = zeros(N_theta,N+1); nu_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1); f_nmt = ones(N_theta,1);
xi_n(:,0+1) = Ar*besselh(r,k*a).*exp(1i*r.*theta);
nu_n(:,0+1) = -Ar*k*a*diff_bessel(2,r,1,k*a).*exp(1i*r.*theta);
f_n = f.*f_n;
xi_n(:,1+1) = Ar*k^1*diff_bessel(2,r,1,k*a).*f_n.*exp(1i*r.*theta);
nu_n(:,1+1) = -f/a.*nu_n(:,1)...
      -Ar*a*k^(1+1).*diff_bessel(2,r,1+1,k*a).*f_n.*exp(1i*r.*theta)...
      -Ar*(2*f).*k^1.*diff_bessel(2,r,1,k*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar*(f_theta/a).*(1i*r).*besselh(r,k*a)...
      .*f_nmo.*exp(1i*r.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);

  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_n(:,n+1) = Ar*k^n*diff_bessel(2,r,n,k*a).*f_n.*exp(1i*r.*theta);
  nu_n(:,n+1) = -f/a.*nu_n(:,n-1+1)...
      -Ar*a*k^(n+1).*diff_bessel(2,r,n+1,k*a).*f_n.*exp(1i*r.*theta)...
      -Ar*(2*f).*k^n.*diff_bessel(2,r,n,k*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar*(f.^2/a)*k^(n-1).*diff_bessel(2,r,n-1,k*a).*f_nmt.*exp(1i*r.*theta)...
      +Ar*(f_theta/a)*k^(n-1).*(1i*r).*diff_bessel(2,r,n-1,k*a)...
      .*f_nmo.*exp(1i*r.*theta);
end

tic;
anp = field_fe_helmholtz_polar_exterior(xi_n,f,k,a,p,N_theta,N);
Gn_fe = dno_fe_helmholtz_polar_exterior(anp,f,f_theta,k,a,p,N_theta,N);
t_fe = toc;
tic;
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi_n,f,f_theta,k,a,b,p,N_theta,N,N_r);
Gn_tfe = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k,a,b,p,N_theta,N,N_r);
t_tfe = toc;

fprintf('  t_fe = %g  t_tfe = %g\n',t_fe,t_tfe);
[relerr,nplot] = compute_errors_2d_polar(nu,Gn_fe,Gn_tfe,Eps,N,N_theta);
% [relerr,nplot] = compute_errors_2d_polar(nu,Gn_fe,Eps,N,N_theta);
make_plots_polar(SavePlots,nplot,relerr);