% test_helmholtz_polar_interior.m
%
% Script to test Helmholtz DNO solvers in polar (interior)
%
% XT 12/16

clear all;
close all;

SavePlots = 0;
RunNumber = 1;
L = 2*pi;
k=2;


if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;
  a = 3; 
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  a = 2.0;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  N_theta = 64;
  a = 2.0;
  N = 16;
end


fprintf('test_helmholtz_polar_interior\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k = %g\n',k);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;

Ar = 1; r = 2; % compute a special wavenumber
A=a+Eps.*f;
xi = Ar*besselj(r,k.*A).*exp(1i*r.*theta);
nu =Ar*((k.*A.*(diff_bessel(1,r,1,k.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k.*A)./A ).*exp(1i*r.*theta)); %DNO
xi_n = zeros(N_theta,N); nu_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
xi_n(:,0+1) = Ar*besselj(r,k*a).*exp(1i*r.*theta);
nu_n(:,0+1) = Ar*k*a*diff_bessel(1,r,1,k*a).*exp(1i*r.*theta);
f_n = f.*f_n;
xi_n(:,1+1) = Ar*k^1*diff_bessel(1,r,1,k*a).*f_n.*exp(1i*r.*theta);
nu_n(:,1+1) = -f/a.*nu_n(:,1)...
      +Ar*a*k^(1+1).*diff_bessel(1,r,1+1,k*a).*f_n.*exp(1i*r.*theta)...
      +Ar*(2*f).*k^1.*diff_bessel(1,r,1,k*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar*(f_theta/a).*(1i*r).*besselj(r,k*a)...
      .*f_nmo.*exp(1i*r.*theta);



for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_n(:,n+1) = Ar*k^n*diff_bessel(1,r,n,k*a).*f_n.*exp(1i*r.*theta);
  nu_n(:,n+1) = -f/a.*nu_n(:,n-1+1)...
      +Ar*a*k^(n+1).*diff_bessel(1,r,n+1,k*a).*f_n.*exp(1i*r.*theta)...
      +Ar*(2*f).*k^n.*diff_bessel(1,r,n,k*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar*(f.^2/a)*k^(n-1).*diff_bessel(1,r,n-1,k*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar*(f_theta/a)*k^(n-1).*(1i*r).*diff_bessel(1,r,n-1,k*a)...
      .*f_nmo.*exp(1i*r.*theta);
end

tic;
dnp = field_fe_helmholtz_polar_interior(xi_n,f,k,a,p,N_theta,N);
% dnp1 = field_fe_helmholtz_polar_interior1(xi,f,k,p,N_theta,N,a);
Gn_fe = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k,a,p,N_theta,N);
t_fe = toc;
% tic;
% Gn_oe = dno_oe_helmholtz(xi,f,p,alphap,betap,eep,eem,N_theta,N);
% Gn_oe = Gn_fe;
% t_oe = toc;
% tic;
% un = field_tfe_helmholtz(xi,f,p,alphap,betap,eep,eem,Dy,a,N_theta,Ny,N);
% Gn_tfe = dno_tfe_helmholtz(un,f,p,alphap,betap,eep,eem,Dy,a,N_theta,Ny,N);
% Gn_tfe = Gn_fe;
% t_tfe = toc;
fprintf('  t_fe = %g\n',t_fe);
% fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
% fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
%     t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);
% 
% [relerr,nplot] = compute_errors_2d(nu,Gn_oe,Gn_fe,Gn_tfe,Eps,N,N_theta);
[relerr,nplot] = compute_errors_2d_polar(nu,Gn_fe,Eps,N,N_theta);

% make_plots_polar(SavePlots,nplot,relerr);