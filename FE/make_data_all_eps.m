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
% Eps_max = 0.01*a;
% name = 'expcos_eps1_VACAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps1_VACAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps1_VACAg';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps1_VACAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos_eps10_VACAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2_eps10_VACAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4_eps10_VACAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps10_VACAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.2*a;
% name = 'expcos_eps20_VACAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps20_VACAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps20_VACAg';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps20_VACAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos_eps40_VACAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps40_VACAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 

%% VACUUM&GOLD
OUT = 'VACUUM';
IN = 'GOLD';
% Eps_max = 0.01*a;
% name = 'expcos_eps1_VACAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps1_VACAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps1_VACAu';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps1_VACAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos_eps10_VACAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2_eps10_VACAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4_eps10_VACAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps10_VACAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.2*a;
% name = 'expcos_eps20_VACAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps20_VACAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps20_VACAu';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps20_VACAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos_eps40_VACAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps40_VACAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% 

%% WATER&SILVER
OUT = 'WATER';
IN = 'SILVER';
% Eps_max = 0.01*a;
% name = 'expcos_eps1_WATERAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps1_WATERAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps1_WATERAg';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps1_WATERAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos_eps10_WATERAg';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2_eps10_WATERAg';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4_eps10_WATERAg';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos810_eps_WATERAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.2*a;
% name = 'expcos_eps20_WATERAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps20_WATERAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps20_WATERAg';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps20_WATERAg';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos_eps40_WATERAg';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps40_WATERAg';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

%% WATER&GOLD
OUT = 'WATER';
IN = 'GOLD';
% Eps_max = 0.01*a;
% name = 'expcos_eps1_WATERAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps1_WATERAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps1_WATERAu';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps1_WATERAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

Eps_max = 0.1*a;
name = 'expcos_eps10_WATERAu';
FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos2_eps10_WATERAu';
FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
name = 'cos4_eps10_WATERAu';
FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps10_WATERAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.2*a;
% name = 'expcos_eps20_WATERAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps20_WATERAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos4_eps20_WATERAu';
% FE_app_pade(M,f4,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos8_eps20_WATERAu';
% FE_app_pade(M,f8,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

% Eps_max = 0.4*a;
% name = 'expcos_eps40_WATERAu';
% FE_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% name = 'cos2_eps40_WATERAu';
% FE_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);