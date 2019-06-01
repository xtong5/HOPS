% save errors
%
% Script to save DNO fe errors in polar (two layers)
%
% XT 11/18

clear all;
close all;
warning off;

%% first load parameters
load('DNO_fe_Eps_0.005.mat');
filename = sprintf('errors_fe_DNO_Eps_0.005.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);

save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% second
load('DNO_fe_Eps_0.01.mat');
filename = sprintf('errors_DNO_fe_Eps_0.01.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% third
load('DNO_fe_Eps_0.05.mat');
filename = sprintf('errors_DNO_fe_Eps_0.05.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 4
load('DNO_fe_Eps_0.1.mat');
filename = sprintf('errors_DNO_fe_Eps_0.1.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 5
load('DNO_fe_Eps_0.005_sing12.mat');
filename = sprintf('errors_DNO_fe_Eps_0.005_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 6
load('DNO_fe_Eps_0.01_sing12.mat');
filename = sprintf('errors_DNO_fe_Eps_0.01_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 7
load('DNO_fe_Eps_0.05_sing12.mat');
filename = sprintf('errors_DNO_fe_Eps_0.05_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 8
load('DNO_fe_Eps_0.1_sing12.mat');
filename = sprintf('errors_DNO_fe_Eps_0.1_sing12.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 9
load('DNO_fe_Eps_0.005_sing16.mat');
filename = sprintf('errors_DNO_fe_Eps_0.005_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 10
load('DNO_fe_Eps_0.01_sing16.mat');
filename = sprintf('errors_DNO_fe_Eps_0.01_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 11
load('DNO_fe_Eps_0.05_sing16.mat');
filename = sprintf('errors_DNO_fe_Eps_0.05_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');

%% 12
load('DNO_fe_Eps_0.1_sing16.mat');
filename = sprintf('errors_DNO_fe_Eps_0.1_sing16.mat');
[relerr_taylorU,relerr_padeU,nplotU] = save_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerr_taylorDNOU,relerr_padeDNOU,nplotDNOU] = save_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
[relerr_taylorW,relerr_padeW,nplotW] = save_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
[relerr_taylorDNOW,relerr_padeDNOW,nplotDNOW] = save_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
save(filename,'Eps','relerr_taylorU','relerr_padeU','nplotU',...
    'relerr_taylorDNOU','relerr_padeDNOU','nplotDNOU',...
    'relerr_taylorW','relerr_padeW','nplotW',...
    'relerr_taylorDNOW','relerr_padeDNOW','nplotDNOW');


