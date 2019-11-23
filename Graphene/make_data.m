% MAKE_DATA.m 
%
% compute and save the date to file 
%
% XT 10/19

clear all
warning('off')

L = 2*pi;
a = 1;
b = 10*a;
N = 16;
M = 101;
N_eps = 101;
N_theta = 64;
lambda_low = 34;
lambda_high = 40;
theta = (L/N_theta)*[0:N_theta-1]';

% name = 'expcos_eps20_VACVAC_35to45_gbar0.025';

OUT = 'VACUUM';
% IN = 'SILVER';
IN = 'VACUUM';
Eps_max = 0.2*a;
mu = 0.4;


% f0 = exp(cos(theta)); name = 'expcos';
% f2 = cos(2*theta); name = 'cos2';    
f4 = cos(4*theta); name = 'cos4';
% f8 = cos(8*theta);
filename = sprintf('%s_eps20_VACVAC_%.0fto%.0f_gbar1_BFPV',name,lambda_low,lambda_high);



tic
% FE_app_TM_pade(f4,N_theta,theta,a,b,N,M,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,name);
FE_app_TM_pade_BFPV(f4,N_theta,theta,a,b,N,M,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,filename);
toc