% data_FE_WATERAg_cos4.m
%
% generate DNO data of exterior layer campared to FE 
% 
% exterior = water, interior = silver, lambda = 0.4125, f = cos(4*theta)
%
% XT 11/17

clear all;
close all;
warning off;
SavePlots = 0;

Mode = 2; 

OUT = 'WATER'; IN = 'SILVER';
a = 0.025;
% N = 24;
N = 16;
N_theta = 64;
Eps = a/5;
  
% lambda = 0.415;
lambda = 0.4275;
[n_u,epsilon_u] = ri_perm(lambda,OUT);
[n_w,epsilon_w] = ri_perm(lambda,IN);

L = 2*pi;
k_0 = L/lambda;
k_u = n_u*k_0; 
k_w = n_w*k_0;



fprintf('data_FE_WATERAg_cos4\n');
fprintf('-------------\n');
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f4 = cos(4*theta);
f4_theta = ifft((1i*p).*fft(f4));

zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u*a.*Sin);
zeta_n(:,0+1) = -Exp; 
psi_n(:,0+1) = (1i*k_u)*a*Sin.*Exp;

for n=1:N
    f_n = f4.*f_n/n;
    if n > 1
        f_nmo = f4.*f_nmo/(n-1);
    end
    zeta_n(:,n+1) = -Exp.*(-1i*k_u)^n.*f_n.*Sin.^n;
    psi_n(:,n+1) = (1i*k_u).*Exp.*(a*(-1i*k_u)^n.*f_n.*...
    Sin.^(n+1)+(f4.*Sin-f4_theta.*cos(theta)).*...
    (-1i*k_u)^(n-1).*f_nmo.*Sin.^(n-1));
 end

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end


% Two-layer scattering by DNO

fprintf('\n\nTwo-layer scattering by DNO\n\n');

tic;
U_n = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f4,f4_theta,tau2,...
    p,k_u,k_w,a,N_theta,N);
apn_fe = field_fe_helmholtz_polar_exterior(U_n,f4,k_u,a,p,N_theta,N);
Gn_fe_u = dno_fe_helmholtz_polar_exterior(apn_fe,f4,f4_theta,k_u,a,p,N_theta,N);
% W_n = U_n - zeta_n;
% dpn_fe = field_fe_helmholtz_polar_interior(W_n,f2,k_w,a,p,N_theta,N);
% Gn_fe_w = dno_fe_helmholtz_polar_interior(dpn_fe,f2,f2_theta,k_w,a,p,N_theta,N);
t_fe = toc;

% fprintf('Press key to compute exterior layer errors...\n');
% pause;

% fprintf('  t_tfe = %g\n',t_tfe);
% fprintf('\nEXTERIOR LAYER\n\n');
% [relerrU,nplotU] = compute_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
% [relerrDNOU,nplotDNOU] = compute_errors_2d_polar(nu_u,Gn_tfe_u,Eps,N,N_theta);

% make_plots_polar(SavePlots,nplotU,relerrU);
% make_plots_polar(SavePlots,nplotDNOU,relerrDNOU);
% fprintf('\n');

% fprintf('Press key to compute interior layer errors...\n');
% pause;

% fprintf('\nINTERIOR LAYER\n\n');
% [relerrW,nplotW] = compute_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
% [relerrDNOW,nplotDNOW] = compute_errors_2d_polar(nu_w,Gn_tfe_w,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotW,relerrW);
% make_plots_polar(SavePlots,nplotDNOW,relerrDNOW);
% fprintf('\n');


filename = sprintf('FE_cos4_eps%g_WATERAg.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'Gn_fe_u','OUT','IN')


