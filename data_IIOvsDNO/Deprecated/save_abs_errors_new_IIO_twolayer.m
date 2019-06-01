% save errors
%
% Script to save new IIO errors in polar (two layers)
%
% XT 5/18

clear all;
close all;
warning off;

%% first load parameters
load('IIO_new_Eps_0.005_Nr32.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.005_Nr32.mat');
L = 2*pi;
theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;
AA=a+Eps.*f;

Ar_u = 2; Ar_w = 1; pp = 2; % take a special wavenumber
 
xi_u = Ar_u*besselh(pp,k_u.*AA).*exp(1i*pp.*theta);
nu_u = Ar_u*((-k_u.*AA.*(diff_besselh(pp,1,k_u.*AA))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k_u.*AA)./AA ).*exp(1i*pp.*theta)); 
xi_w = Ar_w*besselj(pp,k_w.*AA).*exp(1i*pp.*theta);
nu_w = Ar_w*((k_w.*AA.*(diff_besselj(pp,1,k_w.*AA))-...
    1i*pp*Eps.*f_theta.*besselj(pp,k_w.*AA)./AA ).*exp(1i*pp.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar_u*besselh(pp,k_u*a).*exp(1i*pp.*theta);
xi_u_n(:,1+1) = Ar_u*k_u^1*diff_besselh(pp,1,k_u*a).*f_n.*exp(1i*pp.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_besselh(pp,1,k_u*a).*exp(1i*pp.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_besselh(pp,1+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_besselh(pp,1,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar_u*(f_theta/a).*(1i*pp).*besselh(pp,k_u*a).*f_nmo.*exp(1i*pp.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(pp,k_w*a).*exp(1i*pp.*theta);
xi_w_n(:,1+1) = Ar_w*k_w^1*diff_besselj(pp,1,k_w*a).*f_n.*exp(1i*pp.*theta);
nu_w_n(:,0+1) = Ar_w*k_w*a*diff_besselj(pp,1,k_w*a).*exp(1i*pp.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_besselj(pp,1+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_besselj(pp,1,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar_w*(f_theta/a).*(1i*pp).*besselj(pp,k_w*a).*f_nmo.*exp(1i*pp.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_besselh(pp,n,k_u*a).*f_n.*exp(1i*pp.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_besselh(pp,n+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_besselh(pp,n,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_besselh(pp,n-1,k_u*a).*f_nmt.*exp(1i*pp.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*pp).*diff_besselh(pp,n-1,k_u*a)...
      .*f_nmo.*exp(1i*pp.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_besselj(pp,n,k_w*a).*f_n.*exp(1i*pp.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_besselj(pp,n+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_besselj(pp,n,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_besselj(pp,n-1,k_w*a).*f_nmt.*exp(1i*pp.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*pp).*diff_besselj(pp,n-1,k_w*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end


Q_u = 1/sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = 1/sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);

save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');
%% second
load('IIO_new_Eps_0.01_Nr32.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.01_Nr32.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);

save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% third
load('IIO_new_Eps_0.05_Nr32.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.05_Nr32.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 4
load('IIO_new_Eps_0.1_Nr32.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.1_Nr32.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 5
load('IIO_new_Eps_0.005_Nr32_sing12.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.005_Nr32_sing12.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 6
load('IIO_new_Eps_0.01_Nr32_sing12.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.01_Nr32_sing12.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 7
load('IIO_new_Eps_0.05_Nr32_sing12.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.05_Nr32_sing12.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 8
load('IIO_new_Eps_0.1_Nr32_sing12.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.1_Nr32_sing12.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 9
load('IIO_new_Eps_0.005_Nr32_sing16.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.005_Nr32_sing16.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 10
load('IIO_new_Eps_0.01_Nr32_sing16.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.01_Nr32_sing16.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 11
load('IIO_new_Eps_0.05_Nr32_sing16.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.05_Nr32_sing16.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');

%% 12
load('IIO_new_Eps_0.1_Nr32_sing16.mat');
filename = sprintf('abs_errors_IIO_new_Eps_0.1_Nr32_sing16.mat');
[err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
    'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
    'err_taylorW','err_padeW','nplotW',...
    'err_taylorIIOW','err_padeIIOW','nplotIIOW');



% % Two-layer scattering by IIO
% % 
% % fprintf('\n\nTwo-layer scattering by IIO\n\n');
% % 
% % [err_taylorU,err_padeU,nplotU] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
% % [err_taylorIIOU,err_padeIIOU,nplotIIOU] = save_abs_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
% % [err_taylorW,err_padeW,nplotW] = save_abs_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
% % [err_taylorIIOW,err_padeIIOW,nplotIIOW] = save_abs_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
% % 
% % 
% % save(filename,'Eps','err_taylorU','err_padeU','nplotU',...
% %     'err_taylorIIOU','err_padeIIOU','nplotIIOU',...
% %     'err_taylorW','err_padeW','nplotW',...
% %     'err_taylorIIOW','err_padeIIOW','nplotIIOW');

