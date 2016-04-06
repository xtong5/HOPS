clear all;
close all;

RunNumber = 1;
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


% fprintf('test_helmholtz\n');
% fprintf('-------------\n');
% fprintf('RunNumber = %d\n',RunNumber);
% fprintf('alpha = %g  beta = %g\n',alpha,beta);
% fprintf('Eps = %g  a = %g\n',Eps,a);
% fprintf('N_theta = %d  Ny = %d  N = %d\n',N_theta,Ny,N);
% fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;

Ar = -3.0; r = 2; % compute a special index 
xi = Ar*besselh(r,k*a)*exp(li*r.*theta);
A=a+Eps.*f;
nu =(-0.5*k.*A*(besselh(r-1,k.*A)-besselh(r+1,k.*A))+...
    li*r*Eps.f_theta*besselh(r,k.*A)./A ).*exp(li*r.*theta); %DNOß

tic;
apn = field_fe_helmholtz(xi,f,p,alphap,betap,eep,eem,N_theta,N);
Gn_fe = dno_fe_helmholtz(apn,f,p,alphap,betap,eep,eem,N_theta,N);
t_fe = toc;
tic;
%Gn_oe = dno_oe_helmholtz(xi,f,p,alphap,betap,eep,eem,N_theta,N);
Gn_oe = Gn_fe;
t_oe = toc;
tic;
%un = field_tfe_helmholtz(xi,f,p,alphap,betap,eep,eem,Dy,a,N_theta,Ny,N);
%Gn_tfe = dno_tfe_helmholtz(un,f,p,alphap,betap,eep,eem,Dy,a,N_theta,Ny,N);
Gn_tfe = Gn_fe;
t_tfe = toc;
fprintf('  t_fe = %g  t_oe = %g  t_tfe = %g\n',t_oe,t_fe,t_tfe);
fprintf('  t_fe/t_oe = %g  t_tfe/t_fe = %g  t_tfe/t_oe = %g\n',...
    t_fe/t_oe,t_tfe/t_fe,t_tfe/t_oe);

[relerr,nplot] = compute_errors_2d(nu,Gn_oe,Gn_fe,Gn_tfe,Eps,N,N_theta);

% make_plots(SavePlots,nplot,relerr);