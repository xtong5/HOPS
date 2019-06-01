% test: compute and save the date to file 
%
% XT 1/19

L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16; 
M = 21; %test
N_eps = 101; %test
theta = (L/N_theta)*[0:N_theta-1]';

warning('off')

f0 = exp(cos(theta));
f2 = cos(2*theta);
f4 = cos(4*theta);
f8 = cos(8*theta);

%% VACCUM&SILVER
OUT = 'VACUUM';
IN = 'SILVER';
Eps_max = 0.1*a;
name = 'expcos_eps10_VACAg_test';
tic
FE_IIO_app_pade(M,f0,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
toc

% name = 'cos2_eps20_VACAg_test';
% tic
% FE_IIO_app_pade(M,f2,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
% toc