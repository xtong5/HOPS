% test_helmholtz_twolayer_3d.m
%
% Script to test Helmholtz DNO solvers (3d).
%
% DPN 6/7/14
% DPN 1/28/15
% DPN 2/15/15
% DPN 2/24/15
% DPN 2/25/15 ["true" FE]
% DPN 6/21/16

clear all;
clf;
SavePlots = 0;

RunNumber = 1;
Mode = 2;

L1 = 2*pi;
L2 = 2*pi;
alpha1 = 0.1;
alpha2 = 0.2;
beta_u = 1.21;
beta_w = 2.23;
if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  %Nx1 = 64;
  Nx1 = 32;
  Nx1 = 16;
  Nx2 = Nx1;
  %Ny = 24;
  Ny = 32;
  Ny = 16;
  a = 0.5;
  b = 0.5;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  Nx1 = 256;
  Nx2 = Nx1;
  Ny = 64;
  a = 2.0;
  b = 2.0;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  Nx1 = 256;
  Nx2 = Nx1;
  Ny = 64;
  a = 2.0;
  b = 2.0;
  N = 16;
end
k_u = sqrt(alpha1^2 + alpha2^2 + beta_u^2);
k_w = sqrt(alpha1^2 + alpha2^2 + beta_w^2);

fprintf('test_helmholtz_twolayer_3d\n');
fprintf('--------------------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('Mode = %d\n',Mode);
fprintf('alpha1 = %g  alpha2 = %g  beta_u = %g  beta_w = %g\n',alpha1,alpha2,beta_u,beta_w);
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('Nx1 = %d  Nx2 = %d  Ny = %d  N = %d\n',Nx1,Nx2,Ny,N);
fprintf('\n');

x1_pre = (L1/Nx1)*[0:Nx1-1]';
x2_pre = (L2/Nx2)*[0:Nx2-1]';
p1_pre = (2*pi/L1)*[0:Nx1/2-1,-Nx1/2:-1]';
p2_pre = (2*pi/L2)*[0:Nx2/2-1,-Nx2/2:-1]';
x1 = zeros(Nx1,Nx2);
x2 = zeros(Nx1,Nx2);
p1 = zeros(Nx1,Nx2);
p2 = zeros(Nx1,Nx2);

for j1=1:Nx1
  for j2=1:Nx2
    x1(j1,j2) = x1_pre(j1);
    x2(j1,j2) = x2_pre(j2);
    p1(j1,j2) = p1_pre(j1);
    p2(j1,j2) = p2_pre(j2);
  end
end
alpha1p = alpha1 + p1;
alpha2p = alpha2 + p2;
beta_up = 0*alpha1p;
beta_wp = 0*alpha1p;
for j1=1:Nx1
  for j2=1:Nx2
    if(alpha1p(j1,j2)^2 + alpha2p(j1,j2)^2 < k_u^2)
      beta_up(j1,j2) = sqrt(k_u^2 - alpha1p(j1,j2)^2 - alpha2p(j1,j2)^2);
    else
      beta_up(j1,j2) = 1i*sqrt(alpha1p(j1,j2)^2 + alpha2p(j1,j2)^2 - k_u^2);
    end
    if(alpha1p(j1,j2)^2 + alpha2p(j1,j2)^2 < k_w^2)
      beta_wp(j1,j2) = sqrt(k_w^2 - alpha1p(j1,j2)^2 - alpha2p(j1,j2)^2);
    else
      beta_wp(j1,j2) = 1i*sqrt(alpha1p(j1,j2)^2 + alpha2p(j1,j2)^2 - k_w^2);
    end
  end
end
[Dy,y] = cheb(Ny);
eem = exp(-1i*alpha1*x1-1i*alpha2*x2);
eep = exp(1i*alpha1*x1+1i*alpha2*x2);

f = exp(cos(x1)+cos(x2));
f_x1 = -sin(x1).*f;
f_x2 = -sin(x2).*f;

A_u_r = -3.0; r1 = 2; r2 = 2;
alpha1p_r = alpha1p(r1+1,r2+1); alpha2p_r = alpha2p(r1+1,r2+1);
beta_up_r = beta_up(r1+1,r2+1); expia = exp(1i*alpha1p_r*x1+1i*alpha2p_r*x2);
xi_u_n = zeros(Nx1,Nx2,N+1); nu_u_n = zeros(Nx1,Nx2,N+1);
f_n = ones(Nx1,Nx2); f_nmo = ones(Nx1,Nx2);
xi_u_n(:,:,0+1) = A_u_r*expia;
nu_u_n(:,:,0+1) = (-1i*beta_up_r)*A_u_r*expia;
for n=1:N
  f_n = f.*f_n/n;
  if(n>1)
    f_nmo = f.*f_nmo/(n-1);
  end
  xi_u_n(:,:,n+1) = A_u_r*expia.*(f_n*(1i*beta_up_r)^n);
  nu_u_n(:,:,n+1) = (-1i*beta_up_r)*A_u_r*expia.*(f_n*(1i*beta_up_r)^n)...
      +(1i*alpha1p_r)*A_u_r*expia.*(f_x1.*f_nmo*(1i*beta_up_r)^(n-1))...
      +(1i*alpha2p_r)*A_u_r*expia.*(f_x2.*f_nmo*(1i*beta_up_r)^(n-1));
end
xi_u = A_u_r*exp(1i*alpha1p_r*x1+1i*alpha2p_r*x2).*exp(1i*beta_up_r*Eps*f);
nu_u = (-1i*beta_up_r+1i*alpha1p_r*Eps*f_x1+1i*alpha2p_r*Eps*f_x2).*xi_u;

A_w_r = 4.0; r1 = 3; r2 = 1;
alpha1p_r = alpha1p(r1+1,r2+1); alpha2p_r = alpha2p(r1+1,r2+1);
beta_wp_r = beta_wp(r1+1,r2+1); expia = exp(1i*alpha1p_r*x1+1i*alpha2p_r*x2);
xi_w_n = zeros(Nx1,Nx2,N+1); nu_w_n = zeros(Nx1,Nx2,N+1);
f_n = ones(Nx1,Nx2); f_nmo = ones(Nx1,Nx2);
xi_w_n(:,:,0+1) = A_w_r*expia;
nu_w_n(:,:,0+1) = (-1i*beta_wp_r)*A_w_r*expia;
for n=1:N
  f_n = f.*f_n/n;
  if(n>1)
    f_nmo = f.*f_nmo/(n-1);
  end
  xi_w_n(:,:,n+1) = A_w_r*expia.*(f_n*(-1i*beta_wp_r)^n);
  nu_w_n(:,:,n+1) = (-1i*beta_wp_r)*A_w_r*expia.*(f_n*(-1i*beta_wp_r)^n)...
      -(1i*alpha1p_r)*A_w_r*expia.*(f_x1.*f_nmo*(-1i*beta_wp_r)^(n-1))...
      -(1i*alpha2p_r)*A_w_r*expia.*(f_x2.*f_nmo*(-1i*beta_wp_r)^(n-1));
end
xi_w = A_w_r*exp(1i*alpha1p_r*x1+1i*alpha2p_r*x2).*exp(-1i*beta_wp_r*Eps*f);
nu_w = (-1i*beta_wp_r-1i*alpha1p_r*Eps*f_x1-1i*alpha2p_r*Eps*f_x2).*xi_w;

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;

tic;
[apn,dpn] = field_fe_helmholtz_twolayer_3d(zeta_n,psi_n,f,tau2,...
    p1,p2,alpha1p,alpha2p,beta_up,beta_wp,eep,eem,Nx1,Nx2,N);
Gn_fe_u = dno_fe_helmholtz_upper_3d(apn,f,p1,p2,alpha1p,alpha2p,beta_up,eep,eem,Nx1,Nx2,N);
Gn_fe_w = dno_fe_helmholtz_lower_3d(dpn,f,p1,p2,alpha1p,alpha2p,beta_wp,eep,eem,Nx1,Nx2,N);
t_fe = toc;

tic;
[un_u,un_w] = field_tfe_helmholtz_twolayer_3d(zeta_n,psi_n,f,tau2,...
    p1,p2,alpha1p,alpha2p,beta_up,beta_wp,eep,eem,Dy,a,b,Nx1,Nx2,Ny,N);
Gn_tfe_u = dno_tfe_helmholtz_upper_3d(un_u,f,p1,p2,alpha1p,alpha2p,beta_up,...
    eep,eem,Dy,a,Nx1,Nx2,Ny,N);
Gn_tfe_w = dno_tfe_helmholtz_lower_3d(un_w,f,p1,p2,alpha1p,alpha2p,beta_wp,...
    eep,eem,Dy,b,Nx1,Nx2,Ny,N);
t_tfe = toc;

t_oe = t_fe;
fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
    t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);

Gn_oe_u = Gn_fe_u;
[relerr,nplot] = compute_errors_3d(nu_u,Gn_oe_u,Gn_fe_u,Gn_tfe_u,Eps,N,Nx1,Nx2);

fprintf('\n');
fprintf('\nUPPER LAYER\n\n');
fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
fprintf('---------------------------------------------\n');
for n=0:N
  fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
      relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
      relerr(n+1,5),relerr(n+1,6));
end

%make_plots(SavePlots,nplot,relerr);

fprintf('Press key to compute lower layer errors...\n');
pause;

Gn_oe_w = Gn_fe_w;
%Gn_tfe_w = Gn_fe_w;
[relerr,nplot] = compute_errors_3d(nu_w,Gn_oe_w,Gn_fe_w,Gn_tfe_w,Eps,N,Nx1,Nx2);

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
U_n = twolayer_dno_fe_helmholtz_3d(zeta_n,psi_n,f,tau2,...
    p1,p2,alpha1p,alpha2p,beta_up,beta_wp,eep,eem,Nx1,Nx2,N);
apn_fe = field_fe_helmholtz_upper_3d(U_n,f,...
    p1,p2,alpha1p,alpha2p,beta_up,eep,eem,Nx1,Nx2,N);
Gn_fe_u = dno_fe_helmholtz_upper_3d(apn_fe,f,...
    p1,p2,alpha1p,alpha2p,beta_up,eep,eem,Nx1,Nx2,N);
W_n = U_n - zeta_n;
dpn_fe = field_fe_helmholtz_lower_3d(W_n,f,...
    p1,p2,alpha1p,alpha2p,beta_wp,eep,eem,Nx1,Nx2,N);
Gn_fe_w = dno_fe_helmholtz_lower_3d(dpn_fe,f,...
    p1,p2,alpha1p,alpha2p,beta_wp,eep,eem,Nx1,Nx2,N);
t_fe = toc;

tic;
U_n = twolayer_dno_tfe_helmholtz_3d(zeta_n,psi_n,f,tau2,...
    p1,p2,alpha1p,alpha2p,beta_up,beta_wp,eep,eem,Dy,a,b,Nx1,Nx2,Ny,N);
un = field_tfe_helmholtz_upper_3d(U_n,f,...
    p1,p2,alpha1p,alpha2p,beta_up,eep,eem,Dy,a,Nx1,Nx2,Ny,N);
Gn_tfe_u = dno_tfe_helmholtz_upper_3d(un,f,...
    p1,p2,alpha1p,alpha2p,beta_up,eep,eem,Dy,a,Nx1,Nx2,Ny,N);
W_n = U_n - zeta_n;
wn = field_tfe_helmholtz_lower_3d(W_n,f,...
    p1,p2,alpha1p,alpha2p,beta_wp,eep,eem,Dy,b,Nx1,Nx2,Ny,N);
Gn_tfe_w = dno_tfe_helmholtz_lower_3d(wn,f,...
    p1,p2,alpha1p,alpha2p,beta_wp,eep,eem,Dy,b,Nx1,Nx2,Ny,N);
t_tfe = toc;

t_oe = t_fe;
fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
    t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);

Gn_oe_u = Gn_fe_u;
[relerr,nplot] = compute_errors_3d(nu_u,Gn_oe_u,Gn_fe_u,Gn_tfe_u,Eps,N,Nx1,Nx2);

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
[relerr,nplot] = compute_errors_3d(nu_w,Gn_oe_w,Gn_fe_w,Gn_tfe_w,Eps,N,Nx1,Nx2);

fprintf('\n');
fprintf('\nLOWER LAYER\n\n');
fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
fprintf('---------------------------------------------\n');
for n=0:N
  fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
      relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
      relerr(n+1,5),relerr(n+1,6));
end