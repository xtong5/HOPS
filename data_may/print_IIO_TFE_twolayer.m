% print test_IIO_TFE_twolayer.m
%
% Script to test IIO S0=0 and Q0=0 in polar (two layers)
%
% XT 4/18

% clear all;
% close all;
warning off;
SavePlots = 0;

% load('IIO_Eps_0.005_Nr32_eta3.4.mat');
% load('IIO_Eps_0.01_Nr32_eta3.4.mat');
% load('IIO_Eps_0.05_Nr32_eta3.4.mat');
load('IIO_Eps_0.1_Nr32_eta3.4.mat');
% load('IIO_Eps_0.005_Nr32_eta3.4_sing12.mat');
% load('IIO_Eps_0.01_Nr32_eta3.4_sing12.mat');
% load('IIO_Eps_0.05_Nr32_eta3.4_sing12.mat');
% load('IIO_Eps_0.1_Nr32_eta3.4_sing12.mat');
% load('IIO_Eps_0.005_Nr32_eta3.4_sing16.mat');
% load('IIO_Eps_0.01_Nr32_eta3.4_sing16.mat');
% load('IIO_Eps_0.05_Nr32_eta3.4_sing16.mat');
% load('IIO_Eps_0.1_Nr32_eta3.4_sing16.mat');

fprintf('test_IIO_TFE_twolayer\n');
fprintf('-------------\n');
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g  b = %g  c = %g\n',Eps,a,b,c);
fprintf('N_theta = %d N = %d  N_r = %d\n',N_theta,N,N_r);
fprintf('\n');

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


Q_u = 1/sigma_u*nu_u-1i*eta*xi_u;
S_w = 1/sigma_w*nu_w-1i*eta*xi_w;

% Two-layer scattering by IIO

fprintf('\n\nTwo-layer scattering by IIO\n\n');


fprintf('Press key to compute exterior layer errors...\n');
% pause;

% fprintf('  t_tfe = %g\n',t_tfe);
% fprintf('\nEXTERIOR LAYER\n\n');
[relerrU,nplotU] = compute_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerrIIOU,nplotIIOU] = compute_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
% [errIIOU,nplotIIOU] = compute_abserrors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotU,relerrU);
% make_plots_polar(SavePlots,nplotIIOU,relerrIIOU);
% fprintf('\n');

% fprintf('Press key to compute interior layer errors...\n');
% pause;

% fprintf('\nINTERIOR LAYER\n\n');
% [relerrW,nplotW] = compute_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
% [relerrIIOW,nplotIIOW] = compute_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
% [errIIOW,nplotIIOW] = compute_abserrors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotW,relerrW);
% make_plots_polar(SavePlots,nplotIIOW,relerrIIOW);

