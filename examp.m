%p=[1;2;3];n=2;z=1;

clear all;
close all;

L = 2*pi;
k=1;
Eps = 0.02;
N_theta = 64;
a = 0.1; 
N = 16;
 
theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;

Ar = -3.0; r = 2; % compute a special index 
xi = Ar*besselh(r,k*a)*exp(1i*r.*theta);
A=a+Eps.*f;
nu =(-0.5*k.*A*(diff_besselh(r,1,k.*A))+...
    1i*r*Eps.f_theta*besselh(r,k.*A)./A ).*exp(1i*r.*theta); %DNO


