clear all;
close all;

k=2;
Eps = 0.02;
N_theta = 4;
N_r = 2;
a = 1;
b = 2;
N = 16; n=0;
L = 2*pi;

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';
c = (b+a)/(b-a);
w = k*(b-a)/2;
T_p = k* diff_besselh(p,1,b)/besselh(p,b);


%% pre allocation
F_0 = ones(N_r+1,1);
u_n = zeros(N_theta*(N_r+1),1);

f = exp(cos(theta));
% f = cos(theta);
f_theta = ifft( (1i*p).*fft(f) );

[D,x] = cheb(N_r);

r = ((b-a)*x+b+a)/2;

X = b-r;
D_r = (b-a)/2*D;
Dr = diag(r)*D_r;

Diag = diag(x+c);

DD_r = kron(D_r,eye(N_theta));
DDr = kron(Dr,eye(N_theta));

A = diag(kron(X,f));
B = diag(kron(X,f_theta));

D_theta = 1; %?? in fourier space 




AA = [1,2,3]';
BB = [5,4]';
CC=kron(AA,BB);