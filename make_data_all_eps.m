% compute and save the date to file
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
Eps_max = 0.01*a; %0.1*a
M = 11;
N_eps = 2;

warning('off');

theta = (L/N_theta)*[0:N_theta-1]';

f = exp(cos(theta));
name = 'expcos100_eps';
FE_lambda_eps(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% name = 'expcos10_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,10*Eps_max,N_eps,theta,name);
% name = 'expcos5_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,20*Eps_max,N_eps,theta,name);
% name = 'expcos100pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% name = 'expcos10pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,10*Eps_max,N_eps,theta,name);
% name = 'expcos5pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,20*Eps_max,N_eps,theta,name);
% 
% 
% f = cos(2*theta);
% name = 'cos2100_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% name = 'cos210_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,10*Eps_max,N_eps,theta,name);
% name = 'cos25_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,20*Eps_max,N_eps,theta,name);
% name = 'cos2100pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% name = 'cos210pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,10*Eps_max,N_eps,theta,name);
% name = 'cos25pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,20*Eps_max,N_eps,theta,name);
% 
% f = cos(4*theta);
% name = 'cos4100_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% name = 'cos4100pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% 
% f = cos(8*theta);
% name = 'cos8100_eps';
% FE_lambda_eps(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);
% name = 'cos8100pade_eps';
% FE_lambda_eps_pade(M,f,N_theta,a,b,N,Eps_max,N_eps,theta,name);



