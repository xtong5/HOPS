clear all;
close all;

L = 2*pi;
K=1;


if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;

  a = 0.1;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  Nx = 256;
  Ny = 64;
  a = 2.0;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  Nx = 256;
  Ny = 64;
  a = 2.0;
  N = 16;
end
k = sqrt(alpha^2 + beta^2);

fprintf('test_helmholtz\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('alpha = %g  beta = %g\n',alpha,beta);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('Nx = %d  Ny = %d  N = %d\n',Nx,Ny,N);
fprintf('\n');

theta = (L/Nx)*[0:Nx-1]';
p = (2*pi/L)*[0:Nx/2-1,-Nx/2:-1]';
alphap = alpha + p;
betap = 0*alphap;
for j=1:Nx
  if(alphap(j)^2 < k^2)
    betap(j) = sqrt(k^2 - alphap(j)^2);
  else
    betap(j) = 1i*sqrt(alphap(j)^2 - k^2);
  end
end
%[Dy,y] = cheb(Ny);
eem = exp(-1i*alpha*theta);
eep = exp(1i*alpha*theta);

f = exp(cos(theta));
f_x = -sin(theta).*f;

Ar = -3.0; r = 2; alphap_r = alphap(r+1); betap_r = betap(r+1); %???r
xi = Ar*exp(1i*alphap_r*theta).*exp(1i*betap_r*Eps*f); 
nu = (-1i*betap_r+1i*alphap_r*Eps*f_x).*xi; %DNO

tic;
apn = field_fe_helmholtz(xi,f,p,alphap,betap,eep,eem,Nx,N);
Gn_fe = dno_fe_helmholtz(apn,f,p,alphap,betap,eep,eem,Nx,N);
t_fe = toc;
tic;
%Gn_oe = dno_oe_helmholtz(xi,f,p,alphap,betap,eep,eem,Nx,N);
Gn_oe = Gn_fe;
t_oe = toc;
tic;
%un = field_tfe_helmholtz(xi,f,p,alphap,betap,eep,eem,Dy,a,Nx,Ny,N);
%Gn_tfe = dno_tfe_helmholtz(un,f,p,alphap,betap,eep,eem,Dy,a,Nx,Ny,N);
Gn_tfe = Gn_fe;
t_tfe = toc;
fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
    t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);

[relerr,nplot] = compute_errors_2d(nu,Gn_oe,Gn_fe,Gn_tfe,Eps,N,Nx);

% make_plots(SavePlots,nplot,relerr);