% print test_IIO_FE_twolayer.m
%
% Script to test IIO  in polar (two layers)
%
% XT 11/18

% clear all;
% close all;
warning off;
SavePlots = 0;

% load('IIO_fe_Eps_0.005_Nr32_eta3.4.mat');
% load('IIO_fe_Eps_0.01_Nr32_eta3.4.mat');
% load('IIO_fe_Eps_0.05_Nr32_eta3.4.mat');
% load('IIO_fe_Eps_0.1_eta3.4.mat');
% load('IIO_fe_Eps_0.005_Nr32_eta3.4_sing12.mat');
% load('IIO_fe_Eps_0.01_Nr32_eta3.4_sing12.mat');
% load('IIO_fe_Eps_0.05_Nr32_eta3.4_sing12.mat');
% load('IIO_fe_Eps_0.1_Nr32_eta3.4_sing12.mat');
% load('IIO_fe_Eps_0.005_Nr32_eta3.4_sing16.mat');
% load('IIO_fe_Eps_0.01_Nr32_eta3.4_sing16.mat');
% load('IIO_fe_Eps_0.05_Nr32_eta3.4_sing16.mat');
% load('IIO_fe_Eps_0.1_Nr32_eta3.4_sing16.mat');

% load('IIO_new_fe_Eps_0.005.mat');
% load('IIO_new_fe_Eps_0.01.mat');
% load('IIO_new_fe_Eps_0.05.mat');
% load('IIO_new_fe_Eps_0.1.mat');
% load('IIO_new_fe_Eps_0.005_sing12.mat');
% load('IIO_new_fe_Eps_0.01_sing12.mat');
% load('IIO_new_fe_Eps_0.05_sing12.mat');
% load('IIO_new_fe_Eps_0.1_sing12.mat');
load('IIO_new_fe_Eps_0.005_sing16.mat');
% load('IIO_new_fe_Eps_0.01_sing16.mat');
% load('IIO_new_fe_Eps_0.05_sing16.mat');
% load('IIO_new_fe_Eps_0.1_sing16.mat');

fprintf('test_IIO_FE_twolayer\n');
fprintf('-------------\n');
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d \n',N_theta,N);
fprintf('\n');

% Two-layer scattering by IIO

fprintf('\n\nTwo-layer scattering by IIO\n\n');


fprintf('Press key to compute exterior layer errors...\n');
% pause;

% fprintf('  t_fe = %g\n',t_fe);
% fprintf('\nEXTERIOR LAYER\n\n');
[relerrU,nplotU] = compute_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerrIIOU,nplotIIOU] = compute_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
% [errIIOU,nplotIIOU] = compute_abserrors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
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

