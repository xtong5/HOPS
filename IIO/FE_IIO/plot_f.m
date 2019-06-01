% compute and save the date to file
clear all 
close all
a = 0.25;
gbar = 1;
Eps = 2.0/5.0;

L = 2*pi;
N_theta = 128;
theta = (L/N_theta)*[0:N_theta-1]';

mode = 8;

if mode == 1
    f = exp(cos(theta));
end
if mode == 2
    f = cos(2*theta);
end
if mode == 4
    f = cos(4*theta);
end
if mode == 8
    f = cos(8*theta);
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
end

polarplot(theta,gbar+Eps*f);







