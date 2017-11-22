% data_TFE_WATERAg_cos4.m
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
a = 0.025; b = 10*a; c = 0.1*a;
N = 16;
% N = 24;
N_r = 64;
N_theta = 128;
Eps = a/5;
  
% lambda = 0.415;
lambda = 0.4275;
[n_u,epsilon_u] = ri_perm(lambda,OUT);
[n_w,epsilon_w] = ri_perm(lambda,IN);

L = 2*pi;
k_0 = L/lambda;
k_u = n_u*k_0; 
k_w = n_w*k_0;



fprintf('data_TFE_WATERAg_cos4\n');
fprintf('-------------\n');
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g  b = %g  c = %g\n',Eps,a,b,c);
fprintf('N_theta = %d N = %d  N_r = %d\n',N_theta,N,N_r);
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
U_n = twolayer_dno_tfe_helmholtz_polar(zeta_n,psi_n,f4,f4_theta,tau2,...
    p,k_u,k_w,a,b,c,N_theta,N,N_r);
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(U_n,f4,f4_theta,k_u,a,b,p,N_theta,N,N_r);
Gn_tfe_u = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f4,f4_theta,k_u,a,b,p,N_theta,N,N_r);
% W_n = U_n - zeta_n;
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_helmholtz_polar_interior(W_n,f2,f2_theta,k_w,a,c,p,N_theta,N,N_r);
% Gn_tfe_w = dno_tfe_helmholtz_polar_interior(Dr_Wn,Dp_Wn,f2,f2_theta,k_w,a,c,p,N_theta,N,N_r);
t_tfe = toc;


filename = sprintf('TFE_cos4_eps%g_Nr%g_Nt%d_WATERAg.mat',Eps,N_r,N_theta);
save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a',...
    'b','c','Gn_tfe_u','OUT','IN')


