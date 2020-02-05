% MAKE_DATA_IIO.m 
%
% compute and save the date to file 
%
% XT 12/19

clear all
warning('off')

L = 2*pi;
N = 16;
N_lambda = 121;
N_eps = 101;
N_theta = 64;
lambda_low = 34;
% lambda_high = 40;
lambda_high = 39;
theta = (L/N_theta)*[0:N_theta-1]';
Y_p = 1i*3.4.*ones(N_theta,1); Z_p = -1i*3.4.*ones(N_theta,1);

OUT = 'VACUUM';
% IN = 'SILVER';
IN = 'VACUUM';
mu = 0.4;


% f0 = exp(cos(theta)); name = 'expcos';
% f2 = cos(2*theta); name = 'cos2';    
f4 = cos(4*theta); name = 'cos4';
% f8 = cos(8*theta);
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar0025_ALMA',name,lambda_low,lambda_high);
g_bar = 1;  
b = 10*g_bar; Eps_max = 0.2*g_bar;
filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar1_BFPV',name,lambda_low,lambda_high);


tic
% FE_app_IIO_TM_pade_ALMA(f4,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);
FE_app_IIO_TM_pade_BFPV(f4,N_theta,theta,...
    g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);
toc


% g_bar = 0.025; b = 10*g_bar; Eps_max = 0.2*g_bar;
% f2 = cos(2*theta); name = 'cos2';    
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar0025_BFPV',name,lambda_low,lambda_high);
% FE_app_IIO_TM_pade_BFPV(f2,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);
% 

% g_bar = 0.1; b = 10*g_bar; Eps_max = 0.2*g_bar;
% f2 = cos(2*theta); name = 'cos2';    
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar01_BFPV',name,lambda_low,lambda_high);
% FE_app_IIO_TM_pade_BFPV(f2,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);

% g_bar = 1; b = 10*g_bar; Eps_max = 0.2*g_bar;
% f2 = cos(2*theta); name = 'cos2';    
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar1_BFPV',name,lambda_low,lambda_high);
% FE_app_IIO_TM_pade_BFPV(f2,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);

% g_bar = 0.1; b = 10*g_bar; Eps_max = 0.2*g_bar;
% f4 = cos(4*theta); name = 'cos4';    
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar01_BFPV',name,lambda_low,lambda_high);
% FE_app_IIO_TM_pade_BFPV(f4,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);
% 
% g_bar = 0.025; b = 10*g_bar; Eps_max = 0.2*g_bar;
% f4 = cos(4*theta); name = 'cos4';    
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar0025_BFPV',name,lambda_low,lambda_high);
% FE_app_IIO_TM_pade_BFPV(f4,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);
% 
% g_bar = 1; b = 10*g_bar; Eps_max = 0.2*g_bar;
% f4 = cos(4*theta); name = 'cos4';    
% filename = sprintf('IIO_%s_eps20_VACVAC_%.0fto%.0f_gbar1_BFPV',name,lambda_low,lambda_high);
% FE_app_IIO_TM_pade_BFPV(f4,N_theta,theta,...
%     g_bar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda,filename);


