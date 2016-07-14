% test_helmholtz_twolayer.m
%
% Script to test Helmholtz DNO solvers (2d).
%
% DPN 6/7/14
% DPN 1/28/15
% DPN 2/15/15
% DPN 2/23/15
% DPN 2/25/15 ["true" FE]
% DPN 6/21/16

clear all;
clf;
SavePlots = 0;

RunNumber = 1;
Mode = 2;

L = 2*pi;
alpha = 0.1;
beta_u = 1.21;
beta_w = 2.23;
if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  Nx = 32;
  Ny = 16;
  a = 0.5;
  b = 0.5;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  Nx = 256;
  Ny = 64;
  a = 2.0;
  b = 2.0;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  Nx = 256;
  Ny = 64;
  a = 2.0;
  b = 2.0;
  N = 16;
end
k_u = sqrt(alpha^2 + beta_u^2);
k_w = sqrt(alpha^2 + beta_w^2);

fprintf('test_helmholtz_twolayer\n');
fprintf('-----------------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('Mode = %d\n',Mode);
fprintf('alpha = %g  beta_u = %g  beta_w = %g\n',alpha,beta_u,beta_w);
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('Nx = %d  Ny = %d  N = %d\n',Nx,Ny,N);
fprintf('\n');

x = (L/Nx)*[0:Nx-1]';
p = (2*pi/L)*[0:Nx/2-1,-Nx/2:-1]';
alphap = alpha + p;
beta_up = 0*alphap;
beta_wp = 0*alphap;
for j=1:Nx
  if(alphap(j)^2 < k_u^2)
    beta_up(j) = sqrt(k_u^2 - alphap(j)^2);
  else
    beta_up(j) = 1i*sqrt(alphap(j)^2 - k_u^2);
  end
  if(alphap(j)^2 < k_w^2)
    beta_wp(j) = sqrt(k_w^2 - alphap(j)^2);
  else
    beta_wp(j) = 1i*sqrt(alphap(j)^2 - k_w^2);
  end
end
[Dy,y] = cheb(Ny);
eem = exp(-1i*alpha*x);
eep = exp(1i*alpha*x);

f = exp(cos(x));
f_x = -sin(x).*f;

A_u_r = -3.0; r = 2; alphap_r = alphap(r+1); beta_up_r = beta_up(r+1);
xi_u_n = zeros(Nx,N+1); nu_u_n = zeros(Nx,N+1);
f_n = ones(Nx,1); f_nmo = ones(Nx,1);
xi_u_n(:,0+1) = A_u_r*exp(1i*alphap_r*x);
nu_u_n(:,0+1) = (-1i*beta_up_r)*A_u_r*exp(1i*alphap_r*x);
for n=1:N
  f_n = f.*f_n/n;
  if(n>1)
    f_nmo = f.*f_nmo/(n-1);
  end
  xi_u_n(:,n+1) = A_u_r*exp(1i*alphap_r*x).*(f_n*(1i*beta_up_r)^n);
  nu_u_n(:,n+1) = (-1i*beta_up_r)*A_u_r*exp(1i*alphap_r*x).*(f_n*(1i*beta_up_r)^n)...
      +(1i*alphap_r)*A_u_r*exp(1i*alphap_r*x).*(f_x.*f_nmo*(1i*beta_up_r)^(n-1));
end
xi_u = A_u_r*exp(1i*alphap_r*x).*exp(1i*beta_up_r*Eps*f);
nu_u = (-1i*beta_up_r+1i*alphap_r*Eps*f_x).*xi_u;

A_w_r = 4.0; r = 2; alphap_r = alphap(r+1); beta_wp_r = beta_wp(r+1);
xi_w_n = zeros(Nx,N+1); nu_w_n = zeros(Nx,N+1);
f_n = ones(Nx,1); f_nmo = ones(Nx,1);
xi_w_n(:,0+1) = A_w_r*exp(1i*alphap_r*x);
nu_w_n(:,0+1) = (-1i*beta_wp_r)*A_w_r*exp(1i*alphap_r*x);
for n=1:N
  f_n = f.*f_n/n;
  if(n>1)
    f_nmo = f.*f_nmo/(n-1);
  end
  xi_w_n(:,n+1) = A_w_r*exp(1i*alphap_r*x).*(f_n*(-1i*beta_wp_r)^n);
  nu_w_n(:,n+1) = (-1i*beta_wp_r)*A_w_r*exp(1i*alphap_r*x).*(f_n*(-1i*beta_wp_r)^n)...
      -(1i*alphap_r)*A_w_r*exp(1i*alphap_r*x).*(f_x.*f_nmo*(-1i*beta_wp_r)^(n-1));
end
xi_w = A_w_r*exp(1i*alphap_r*x).*exp(-1i*beta_wp_r*Eps*f);
nu_w = (-1i*beta_wp_r-1i*alphap_r*Eps*f_x).*xi_w;

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;

tic;
[apn,dpn] = field_fe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Nx,N);
Gn_fe_u = dno_fe_helmholtz_upper(apn,f,p,alphap,beta_up,eep,eem,Nx,N);
Gn_fe_w = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,N);
t_fe = toc;

tic;
[un_u,un_w] = field_tfe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Dy,a,b,Nx,Ny,N);
Gn_tfe_u = dno_tfe_helmholtz_upper(un_u,f,p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,N);
Gn_tfe_w = dno_tfe_helmholtz_lower(un_w,f,p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,N);
t_tfe = toc;

t_oe = t_fe;
fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
    t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);

Gn_oe_u = Gn_fe_u;
[relerr,nplot] = compute_errors_2d(nu_u,Gn_oe_u,Gn_fe_u,Gn_tfe_u,Eps,N,Nx);

fprintf('\n');
fprintf('\nUPPER LAYER\n\n');
fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
fprintf('---------------------------------------------\n');
for n=0:N
  fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
      relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
      relerr(n+1,5),relerr(n+1,6));
end

fprintf('Press key to compute lower layer errors...\n');
pause;

Gn_oe_w = Gn_fe_w;
[relerr,nplot] = compute_errors_2d(nu_w,Gn_oe_w,Gn_fe_w,Gn_tfe_w,Eps,N,Nx);

fprintf('\n');
fprintf('\nLOWER LAYER\n\n');
fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
fprintf('---------------------------------------------\n');
for n=0:N
  fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
      relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
      relerr(n+1,5),relerr(n+1,6));
end

%make_plots(SavePlots,nplot,relerr);

%
% Two-layer scattering by DNO
%

fprintf('Press key to compute by DNO...\n');
pause;

fprintf('\n\nTwo-layer scattering by DNO\n\n');

tic;
U_n = dno_fe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Nx,N);
apn_fe = field_fe_helmholtz_upper(U_n,f,p,alphap,beta_up,eep,eem,Nx,N);
Gn_fe_u = dno_fe_helmholtz_upper(apn_fe,f,p,alphap,beta_up,eep,eem,Nx,N);
W_n = U_n - zeta_n;
dpn_fe = field_fe_helmholtz_lower(W_n,f,p,alphap,beta_wp,eep,eem,Nx,N);
Gn_fe_w = dno_fe_helmholtz_lower(dpn_fe,f,p,alphap,beta_wp,eep,eem,Nx,N);
t_fe = toc;

tic;
U_n = dno_tfe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Dy,a,b,Nx,Ny,N);
un = field_tfe_helmholtz_upper(U_n,f,...
    p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,N);
Gn_tfe_u = dno_tfe_helmholtz_upper(un,f,...
    p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,N);
W_n = U_n - zeta_n;
wn = field_tfe_helmholtz_lower(W_n,f,...
    p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,N);
Gn_tfe_w = dno_tfe_helmholtz_lower(wn,f,...
    p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,N);
t_tfe = toc;

t_oe = t_fe;
fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
    t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);

Gn_oe_u = Gn_fe_u;
[relerr,nplot] = compute_errors_2d(nu_u,Gn_oe_u,Gn_fe_u,Gn_tfe_u,Eps,N,Nx);

fprintf('\n');
fprintf('\nUPPER LAYER\n\n');
fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
fprintf('---------------------------------------------\n');
for n=0:N
  fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
      relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
      relerr(n+1,5),relerr(n+1,6));
end

fprintf('Press key to compute lower layer errors...\n');
pause;

Gn_oe_w = Gn_fe_w;
[relerr,nplot] = compute_errors_2d(nu_w,Gn_oe_w,Gn_fe_w,Gn_tfe_w,Eps,N,Nx);

fprintf('\n');
fprintf('\nLOWER LAYER\n\n');
fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
fprintf('---------------------------------------------\n');
for n=0:N
  fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
      relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
      relerr(n+1,5),relerr(n+1,6));
end