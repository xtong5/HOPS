% compute and save the date to file 
L = 2*pi;
a = 0.025;
b = 10*a;
M = 41;
N_eps = 201;
N_theta = 64;
theta = (L/N_theta)*[0:N_theta-1]';

warning('off')

f0 = exp(cos(theta));
% f2 = cos(2*theta);
% f4 = cos(4*theta);
% f8 = cos(8*theta);

% %% N20
% N_theta = 64;N = 20;
% OUT = 'VACUUM';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_VACAgN20';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_WATERAgN20';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);


% %% N24
% N_theta = 64;N = 24;
% OUT = 'VACUUM';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_VACAgN24';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_WATERAgN24';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);


% %% N_theta96 N16
% N_theta = 96;N = 16;
% theta = (L/N_theta)*[0:N_theta-1]';
% OUT = 'VACUUM';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_VACAgNt96N16';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_WATERAgNt96N16';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 
% 
% %% N_theta96 N20
% N_theta = 96;N = 20;
% OUT = 'VACUUM';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_VACAgNt96N20';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_WATERAgNt96N20';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 
% 
% %% N_theta96 N24
% N_theta = 96;N = 24;
% OUT = 'VACUUM';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_VACAgNt96N24';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_WATERAgNt96N24';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

%% N16
N = 16;
OUT = 'VACUUM';
IN = 'SILVER';
Eps_max = 0.2*a;
name = 'expcos_eps20_VACAg_test_dno';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.1*a;
% name = 'cos810_eps_WATERAgN16';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 
% N = 16;
% OUT = 'VACUUM';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos820_eps_VACAgN16';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos820_eps_WATERAgN16';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 
% 
% OUT = 'WATER';
% IN = 'SILVER';
% Eps_max = 0.2*a;
% name = 'cos45_eps_WATERAgN16';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
