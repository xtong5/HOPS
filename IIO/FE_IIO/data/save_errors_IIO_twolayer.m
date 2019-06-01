% save errors
%
% Script to save IIO errors in polar (two layers)
%
% XT 11/18

clear all;
close all;
warning off;

%% 1
load('IIO_fe_Eps_0.005_eta3.4.mat');
filename = sprintf('errors_IIO_fe_Eps_0.005.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 2
load('IIO_fe_Eps_0.01_eta3.4.mat');
filename = sprintf('errors_IIO_fe_Eps_0.01.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 3
load('IIO_fe_Eps_0.05_eta3.4.mat');
filename = sprintf('errors_IIO_fe_Eps_0.05.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 4
load('IIO_fe_Eps_0.1_eta3.4.mat');
filename = sprintf('errors_IIO_fe_Eps_0.1.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 5
load('IIO_fe_Eps_0.005_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_fe_Eps_0.005_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 6
load('IIO_fe_Eps_0.01_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_fe_Eps_0.01_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 7
load('IIO_fe_Eps_0.05_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_fe_Eps_0.05_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 8
load('IIO_fe_Eps_0.1_eta3.4_sing12.mat');
filename = sprintf('errors_IIO_fe_Eps_0.1_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 9
load('IIO_fe_Eps_0.005_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_fe_Eps_0.005_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 10
load('IIO_fe_Eps_0.01_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_fe_Eps_0.01_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 11
load('IIO_fe_Eps_0.05_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_fe_Eps_0.05_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 12
load('IIO_fe_Eps_0.1_eta3.4_sing16.mat');
filename = sprintf('errors_IIO_fe_Eps_0.1_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');


%% NEW IIO
%% 1
load('IIO_new_fe_Eps_0.005.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.005.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 2
load('IIO_new_fe_Eps_0.01.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.01.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 3
load('IIO_new_fe_Eps_0.05.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.05.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 4
load('IIO_new_fe_Eps_0.1.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.1.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 5
load('IIO_new_fe_Eps_0.005_sing12.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.005_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 6
load('IIO_new_fe_Eps_0.01_sing12.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.01_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 7
load('IIO_new_fe_Eps_0.05_sing12.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.05_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 8
load('IIO_new_fe_Eps_0.1_sing12.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.1_sing12.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 9
load('IIO_new_fe_Eps_0.005_sing16.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.005_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 10
load('IIO_new_fe_Eps_0.01_sing16.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.01_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 11
load('IIO_new_fe_Eps_0.05_sing16.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.05_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

%% 12
load('IIO_new_fe_Eps_0.1_sing16.mat');
filename = sprintf('errors_IIO_new_fe_Eps_0.1_sing16.mat');
[relerr_taylorIU,relerr_padeIU,nplotIU] = save_errors_2d_polar(I_u,I_u_n,Eps,N,N_theta);
[relerr_taylorIIOU,relerr_padeIIOU,nplotIIOU] = save_errors_2d_polar(Q_u,Q_u_n,Eps,N,N_theta);
[relerr_taylorIW,relerr_padeIW,nplotIW] = save_errors_2d_polar(I_w,I_w_n,Eps,N,N_theta);
[relerr_taylorIIOW,relerr_padeIIOW,nplotIIOW] = save_errors_2d_polar(S_w,S_w_n,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorIU','relerr_padeIU','nplotIU',...
    'relerr_taylorIIOU','relerr_padeIIOU','nplotIIOU',...
    'relerr_taylorIW','relerr_padeIW','nplotIW',...
    'relerr_taylorIIOW','relerr_padeIIOW','nplotIIOW');

