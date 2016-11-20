% compute and save the date to file
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
Eps = 1; 
M = 201;


theta = (L/N_theta)*[0:N_theta-1]';

ff = zeros(N_theta/4,1);
for i=0:N_theta/8
    ff(i+1) = a/cos(theta(i+1));
end
for i=N_theta/8+1:N_theta/8*2-1
    ff(i+1) = a/sin(theta(i+1));
end
f = [ff; ff; ff; ff]-a;
p = [0:N_theta/2-1,-N_theta/2:-1]';
f_theta = ifft(1i*p.*fft(f));
    
% polarplot(theta,f)
name = 'square';


FE_application(M,f,N_theta,a,b,N,Eps,theta,name);

