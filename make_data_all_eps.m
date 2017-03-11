% compute and save the date to file 
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
M = 201;
N_eps = 201;
theta = (L/N_theta)*[0:N_theta-1]';

warning('off')

f0 = exp(cos(theta));
f2 = cos(2*theta);
f4 = cos(4*theta);
% f8 = cos(8*theta);

%% VACCUM&SILVER
OUT = 'VACUUM';
IN = 'SILVER';
Eps_max = 0.01*a;
name = 'expcos100_eps_VACAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2100_eps_VACAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4100_eps_VACAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8100_eps_VACAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos10_eps_VACAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos210_eps_VACAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos410_eps_VACAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos810_eps_VACAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.2*a;
name = 'expcos5_eps_VACAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos25_eps_VACAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos45_eps_VACAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos85_eps_VACAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos52_eps_VACAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos252_eps_VACAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 

%% VACUUM&GOLD
OUT = 'VACUUM';
IN = 'GOLD';
Eps_max = 0.01*a;
name = 'expcos100_eps_VACAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2100_eps_VACAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4100_eps_VACAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8100_eps_VACAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos10_eps_VACAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos210_eps_VACAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos410_eps_VACAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos810_eps_VACAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.2*a;
name = 'expcos5_eps_VACAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos25_eps_VACAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos45_eps_VACAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos85_eps_VACAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos52_eps_VACAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos252_eps_VACAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 

%% WATER&SILVER
OUT = 'WATER';
IN = 'SILVER';
Eps_max = 0.01*a;
name = 'expcos100_eps_WATERAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2100_eps_WATERAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4100_eps_WATERAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8100_eps_WATERAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos10_eps_WATERAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos210_eps_WATERAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos410_eps_WATERAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos810_eps_WATERAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.2*a;
name = 'expcos5_eps_WATERAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos25_eps_WATERAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos45_eps_WATERAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos85_eps_WATERAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos52_eps_WATERAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos252_eps_WATERAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

%% WATER&GOLD
OUT = 'WATER';
IN = 'GOLD';
Eps_max = 0.01*a;
name = 'expcos100_eps_WATERAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2100_eps_WATERAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4100_eps_WATERAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8100_eps_WATERAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos10_eps_WATERAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos210_eps_WATERAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos410_eps_WATERAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos810_eps_WATERAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.2*a;
name = 'expcos5_eps_WATERAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos25_eps_WATERAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos45_eps_WATERAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos85_eps_WATERAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos52_eps_WATERAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos252_eps_WATERAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);