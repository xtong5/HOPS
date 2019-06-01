% compute and save the date to file 
%
% XT 12/18
% XT 2/19

L = 2*pi;
% N_theta = 64;
N_theta = 128;
a = 0.025;
b = 10*a;
N = 16; 
% M = 41; %test
% N_eps = 101; %test
M = 201; 
N_eps = 201; 
theta = (L/N_theta)*[0:N_theta-1]';

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
% name = 'expcos_eps0_VACAg';
% FE_IIO_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps10_VACAg';
% FE_IIO_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps10_VACAg';
% FE_IIO_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.2*a;
% name = 'expcos_eps20_VACAg';
% FE_IIO_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps20_VACAg';
% FE_IIO_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4_eps20_VACAg_Nth128';
FE_IIO_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);


%% WATER&SILVER
OUT = 'WATER';
IN = 'SILVER';

% Eps_max = 0.1*a;
% name = 'expcos_eps10_WATERAg';
% FE_IIO_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps10_WATERAg';
% FE_IIO_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps10_WATERAg';
% FE_IIO_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.2*a;
% name = 'expcos_eps20_WATERAg';
% FE_IIO_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps20_WATERAg';
% FE_IIO_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4_eps20_WATERAg_Nth128';
FE_IIO_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

toc
