%
% XT 11/18

clear all;
close all;
warning off;
SavePlots = 0;
Mode = 2; %check 

L = 2*pi;
lambda = 0.45;
n_u = 1;
n_w = 2.5;
k_0 = L/lambda;
k_u = n_u*k_0; 
k_w = n_w*k_0;
% k_w = 5.13562230184068;
if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = 1/(lambda*k_u/L)^2;
  sigma_w = 1/(lambda*k_w/L)^2;
end

N_theta = 64;
N = 16;
a = 0.025; 
Eps = 0.02;
b = 10*a;


fprintf('test_farfield\n');
fprintf('-------------\n');
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;

Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

compute_init;
U_far = Ar_u*besselh(pp,k_u.*b).*exp(1i*pp.*theta); % true farfield 

% Z_p = -1i*3.4.*ones(N_theta,1);
% Y_p = 1i*3.4.*ones(N_theta,1);

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
% Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
% dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
% S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
Un_far = IIO_fe_farfield(anp,k_u,a,b,p,N_theta,N);
t_fe = toc;

% norm
for j=1:N_theta
    i = floor(N/2);
    B_far(j) = padesum(Un_far(j,:).',Eps,i);
end
B_norm = norm(B_far,2)/sqrt(N_theta);
       
fprintf('  t_fe = %g\n',t_fe);
fprintf('\nFarfield\n\n');
[relerrUfar,nplotUfar] = compute_errors_2d_polar(U_far,Un_far,Eps,N,N_theta);
fprintf('\n');


%% Eps=0
Eps = 0;
compute_init;
tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp0 = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Un_far0 = IIO_fe_farfield(anp0,k_u,a,b,p,N_theta,N);
t_fe0 = toc;
% norm
for j=1:N_theta
    i = floor(N/2);
    B_far0(j) = padesum(Un_far0(j,:).',Eps,i);
end
B_norm0 = norm(B_far0,2)/sqrt(N_theta);

fprintf('  t_fe = %g\n',t_fe0);
fprintf('\nFarfield\n\n');
[relerrUfar0,nplotUfar0] = compute_errors_2d_polar(U_far,Un_far0,Eps,N,N_theta);
fprintf('\n');

fprintf('B_norm = %g   B_norm0 = %g  \n',B_norm,B_norm0);






