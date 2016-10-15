% compute and save the date to file
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
Eps = 0.01*a; %0.1*a
M = 401;

theta = (L/N_theta)*[0:N_theta-1]';

mode = 0; % choose functions

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
if mode == 0
   ff = zeros(N_theta/4,1);
   for i=0:N_theta/8
       ff(i+1) = a/cos(theta(i+1));
   end
   for i=N_theta/8+1:N_theta/8*2-1
       ff(i+1) = a/sin(theta(i+1));
   end
   f = [ff; ff; ff; ff];
   name = 'square';
end

% f = a.*(theta == 0)+ a./cos(theta).*...
%     ( (theta > 0 & theta <= L/8) | (theta > 7*L/4 & theta <= L) )...
%     + a./sin(theta).* ( (theta > L/8 & theta <= 3*L/8) )...
%     - a./cos(theta).* ( (theta > 3*L/8 & theta <= 5*L/8) )...
%     - a./sin(theta).* ( (theta > 5*L/8 & theta <= 7*L/8) ); 

FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

