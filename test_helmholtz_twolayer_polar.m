% clear all;
% clf;
% SavePlots = 0;

RunNumber = 1;
Mode = 2; %check 

L = 2*pi;
k_u = 1; %diff k?????
k_w = 1;

if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;
  a = 1; 
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

fprintf('test_helmholtz_twolayer_polar\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;
A=a+Eps.*f;

Ar_u = 1; r = 2; % take a special wavenumber
Ar_w = 1; %r = 2; % take a special wavenumber ???same

xi_u = Ar_u*besselh(r,k_u.*A).*exp(1i*r.*theta);
nu_u = Ar_u*((-k_u.*A.*(diff_besselh(r,1,k_u.*A))+...
    1i*r*Eps.*f_theta.*besselh(r,k_u.*A)./A ).*exp(1i*r.*theta)); 
xi_w = Ar_w*besselj(r,k_w.*A).*exp(1i*r.*theta);
nu_w = Ar_w*((k_w.*A.*(diff_besselj(r,1,k_w.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k_w.*A)./A ).*exp(1i*r.*theta)); 
xi_u_n = zeros(N_theta,N); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
xi_u_n(:,0+1) = Ar_u*besselh(r,k_u*a).*exp(1i*r.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_besselh(r,1,k_u*a).*exp(1i*r.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(r,k_w*a).*exp(1i*r.*theta);
nu_w_n(:,0+1) = -Ar_w*k_w*a*diff_besselh(r,1,k_w*a).*exp(1i*r.*theta);


for n=1:N
  f_n = f.*f_n/n;
  if(n>1)
    f_nmo = f.*f_nmo/(n-1);
  end
  if(n>2)
    f_nmo = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_besselh(r,n,k_u*a).*f_n.*exp(1i*r.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_besselh(r,n+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_besselh(r,n,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_besselh(r,n-1,k_u*a).*f_nmt.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*r).*diff_besselh(r,n-1,k_u*a)...
      .*f_nmo.*exp(1i*r.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_besselj(r,n,k_w*a).*f_n.*exp(1i*r.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_besselj(r,n+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_besselj(r,n,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_besselj(r,n-1,k_w*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*r).*diff_besselj(r,n-1,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);
end


if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;
% [anp] = field_fe_helmholtz_twolayer_polar(zeta_n,psi_n,f,f_theta,tau2,...
%     p,k_u,k_w,a,N_theta,N);
[U_n] = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N);


% tic;
% [apn,dpn] = field_fe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
%     p,alphap,beta_up,beta_wp,eep,eem,Nx,N);
% Gn_fe_u = dno_fe_helmholtz_upper(apn,f,p,alphap,beta_up,eep,eem,Nx,N);
% Gn_fe_w = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,N);
% t_fe = toc;
% 
% tic;
% [un_u,un_w] = field_tfe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
%     p,alphap,beta_up,beta_wp,eep,eem,Dy,a,b,Nx,Ny,N);
% Gn_tfe_u = dno_tfe_helmholtz_upper(un_u,f,p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,N);
% Gn_tfe_w = dno_tfe_helmholtz_lower(un_w,f,p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,N);
% t_tfe = toc;
% 
% t_oe = t_fe;
% fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
% fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
%     t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);
% 
% Gn_oe_u = Gn_fe_u;
% [relerr,nplot] = compute_errors_2d(nu_u,Gn_oe_u,Gn_fe_u,Gn_tfe_u,Eps,N,Nx);
% 
% fprintf('\n');
% fprintf('\nUPPER LAYER\n\n');
% fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
% fprintf('---------------------------------------------\n');
% for n=0:N
%   fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
%       relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
%       relerr(n+1,5),relerr(n+1,6));
% end
% 
% fprintf('Press key to compute lower layer errors...\n');
% pause;
% 
% Gn_oe_w = Gn_fe_w;
% [relerr,nplot] = compute_errors_2d(nu_w,Gn_oe_w,Gn_fe_w,Gn_tfe_w,Eps,N,Nx);
% 
% fprintf('\n');
% fprintf('\nLOWER LAYER\n\n');
% fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
% fprintf('---------------------------------------------\n');
% for n=0:N
%   fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
%       relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
%       relerr(n+1,5),relerr(n+1,6));
% end
% 
% %make_plots(SavePlots,nplot,relerr);
% 
% %
% % Two-layer scattering by DNO
% %
% 
% fprintf('Press key to compute by DNO...\n');
% pause;
% 
% fprintf('\n\nTwo-layer scattering by DNO\n\n');
% 
% tic;
% U_n = dno_fe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
%     p,alphap,beta_up,beta_wp,eep,eem,Nx,N);
% apn_fe = field_fe_helmholtz_upper(U_n,f,p,alphap,beta_up,eep,eem,Nx,N);
% Gn_fe_u = dno_fe_helmholtz_upper(apn_fe,f,p,alphap,beta_up,eep,eem,Nx,N);
% W_n = U_n - zeta_n;
% dpn_fe = field_fe_helmholtz_lower(W_n,f,p,alphap,beta_wp,eep,eem,Nx,N);
% Gn_fe_w = dno_fe_helmholtz_lower(dpn_fe,f,p,alphap,beta_wp,eep,eem,Nx,N);
% t_fe = toc;
% 
% tic;
% U_n = dno_tfe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
%     p,alphap,beta_up,beta_wp,eep,eem,Dy,a,b,Nx,Ny,N);
% un = field_tfe_helmholtz_upper(U_n,f,...
%     p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,N);
% Gn_tfe_u = dno_tfe_helmholtz_upper(un,f,...
%     p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,N);
% W_n = U_n - zeta_n;
% wn = field_tfe_helmholtz_lower(W_n,f,...
%     p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,N);
% Gn_tfe_w = dno_tfe_helmholtz_lower(wn,f,...
%     p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,N);
% t_tfe = toc;
% 
% t_oe = t_fe;
% fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
% fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
%     t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);
% 
% Gn_oe_u = Gn_fe_u;
% [relerr,nplot] = compute_errors_2d(nu_u,Gn_oe_u,Gn_fe_u,Gn_tfe_u,Eps,N,Nx);
% 
% fprintf('\n');
% fprintf('\nUPPER LAYER\n\n');
% fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
% fprintf('---------------------------------------------\n');
% for n=0:N
%   fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
%       relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
%       relerr(n+1,5),relerr(n+1,6));
% end
% 
% fprintf('Press key to compute lower layer errors...\n');
% pause;
% 
% Gn_oe_w = Gn_fe_w;
% [relerr,nplot] = compute_errors_2d(nu_w,Gn_oe_w,Gn_fe_w,Gn_tfe_w,Eps,N,Nx);
% 
% fprintf('\n');
% fprintf('\nLOWER LAYER\n\n');
% fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
% fprintf('---------------------------------------------\n');
% for n=0:N
%   fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
%       relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
%       relerr(n+1,5),relerr(n+1,6));
% end