% compute and save the date to file
L = 2*pi;
N_theta = 16;
a = 0.025;
b = 10*a;
N = 16;
Eps = 0.01*a;

theta = (L/N_theta)*[0:N_theta-1]';
% f = exp(cos(theta));
f = cos(4*theta);
% f = cos(2*theta);
% f = cos(8*theta);


FE_application(M,f,N_theta,a,b,N,Eps,L,theta);
