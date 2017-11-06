% compute and save the date to file
% data given by different lambda and functions (fixed epsilon)

L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
% Eps = 0.01*a; %0.1*a
Eps = 0.01*a;
M = 401;

warning('off');

theta = (L/N_theta)*[0:N_theta-1]';

f = exp(cos(theta));
name = 'expcos100pade';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);


f = cos(2*theta);
name = 'cos2100pade';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(4*theta);
name = 'cos4100pade';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(8*theta);
name = 'cos8100pade';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);


