% compute and save the date to file
% data given by different lambda for one function (fixed lambda)

L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
Eps = 0.01*a; %0.1*a
M = 1;

theta = (L/N_theta)*[0:N_theta-1]';

mode = 1; % choose functions

if mode == 1
    f = exp(cos(theta));
    name = 'expcos';
end
if mode == 2
    f = cos(2*theta);
    name = 'cos2';
end
if mode == 4
    f = cos(4*theta);
    name = 'cos4';
end
if mode == 8
    f = cos(8*theta);
    name = 'cos8';
end

FE_applicationTaylor(M,f,N_theta,a,b,N,Eps,theta,name);
%FE_applicationPade(M,f,N_theta,a,b,N,Eps,theta,name);

