% compute and save the date to file 
% 
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
Eps_max = 0.01*a; %0.1*a
M = 201;
N_eps = 101;

theta = (L/N_theta)*[0:N_theta-1]';

mode = 1; % choose functions

if mode == 1
    f = exp(cos(theta));
    name = 'expcos100_eps';
    OUT = 'VACUUM';
    IN = 'SILVER';
end
if mode == 2
    f = cos(2*theta);
    name = 'cos2100_eps';
end
if mode == 4
    f = cos(4*theta);
    name = 'cos4100_eps';
end
if mode == 8
    f = cos(8*theta);
    name = 'cos8100_eps';
end


FE_app_pade(M,f,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);
