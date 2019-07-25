% compute and save the date to file 
%
% XT 6/19

L = 2*pi;
N_theta = 32;
a = 0.025;
b = 10*a;
c = 0.1*a;
N = 8; 
N_r = 16;
% M = 41; %test
% N_eps = 101; %test
M = 201; 
N_eps = 201; 
theta = (L/N_theta)*[0:N_theta-1]';
eta = 3.4;
Y_p = 1i*eta*ones(N_theta,1);
Z_p = -1i*eta*ones(N_theta,1);

warning('off')

f0 = exp(cos(theta));
f2 = cos(2*theta);
f4 = cos(4*theta);
% f8 = cos(8*theta);

tic

%% VACCUM&SILVER
OUT = 'VACUUM';
IN = 'SILVER';


% Eps_max = 0.1*a;
% name = 'expcos_eps10_VACAg';
% TFE_IIO_app_pade(M,f0,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos2_eps10_VACAg';
% TFE_IIO_app_pade(M,f2,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos4_eps10_VACAg';
% TFE_IIO_app_pade(M,f4,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);

Eps_max = 0.2*a;
name = 'expcos_eps20_VACAg';
TFE_IIO_app_pade(M,f0,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
name = 'cos2_eps20_VACAg';
TFE_IIO_app_pade(M,f2,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos4_eps20_VACAg';
% TFE_IIO_app_pade(M,f4,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);


%% WATER&SILVER
% OUT = 'WATER';
% IN = 'SILVER';

% Eps_max = 0.1*a;
% name = 'expcos_eps10_WATERAg';
% TFE_IIO_app_pade(M,f0,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos2_eps10_WATERAg';
% TFE_IIO_app_pade(M,f2,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos4_eps10_WATERAg';
% TFE_IIO_app_pade(M,f4,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);

% Eps_max = 0.2*a;
% name = 'expcos_eps20_WATERAg';
% TFE_IIO_app_pade(M,f0,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos2_eps20_WATERAg';
% TFE_IIO_app_pade(M,f2,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);
% name = 'cos4_eps20_WATERAg_Nth128';
% TFE_IIO_app_pade(M,f4,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,Y_p,Z_p,name);

toc
