% save errors
%
% Script to save IIO errors in polar (two layers)
%
% XT 5/18

clear all;
close all;
warning off;

%% 1
load('IIO_Eps_0.005_Nr32_eta3.4.mat');
filename = sprintf('errors_IIO_Eps_0.005_Nr32.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 2
load('IIO_Eps_0.01_Nr32_eta3.4.mat');
filename = sprintf('errors_IIO_Eps_0.01_Nr32.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 3
load('IIO_Eps_0.05_Nr32_eta3.4.mat');
filename = sprintf('errors_IIO_Eps_0.05_Nr32.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 4
load('IIO_Eps_0.1_Nr32_eta3.4.mat');
filename = sprintf('errors_IIO_Eps_0.1_Nr32.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 5
load('IIO_Eps_0.005_Nr32_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_Eps_0.005_Nr32_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 6
load('IIO_Eps_0.01_Nr32_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_Eps_0.01_Nr32_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 7
load('IIO_Eps_0.05_Nr32_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_Eps_0.05_Nr32_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 8
load('IIO_Eps_0.1_Nr32_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_Eps_0.1_Nr32_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 9
load('IIO_Eps_0.005_Nr32_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_Eps_0.005_Nr32_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 10
load('IIO_Eps_0.01_Nr32_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_Eps_0.01_Nr32_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 11
load('IIO_Eps_0.05_Nr32_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_Eps_0.05_Nr32_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 12
load('IIO_Eps_0.1_Nr32_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_Eps_0.1_Nr32_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,Un,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Qn_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,Wn,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,Sn_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');
