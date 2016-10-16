% compute and save the date to file
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
% Eps = 0.01*a; %0.1*a
Eps = 0.1*a;
M = 401;

theta = (L/N_theta)*[0:N_theta-1]';

f = exp(cos(theta));
name = 'expcos10';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);


f = cos(2*theta);
name = 'cos210';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(4*theta);
name = 'cos410';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(8*theta);
name = 'cos810';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);


ff = zeros(N_theta/4,1);
for i=0:N_theta/8
    ff(i+1) = a/cos(theta(i+1));
end
for i=N_theta/8+1:N_theta/8*2-1
    ff(i+1) = a/sin(theta(i+1));
end
f = [ff; ff; ff; ff];
name = 'square10';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);


