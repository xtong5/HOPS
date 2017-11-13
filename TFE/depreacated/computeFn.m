clear all;
close all;

k=2;
Eps = 0.02;
N_theta = 64;
N_r = 16;
a = 1;
b = 2;
d = b-a;
N = 16; 
L = 2*pi;

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

c = (b+a)/d;
w = k*d/2;
T_p = k* diff_besselh(p,1,b)/besselh(p,b);


%% pre allocation
Un = zeros((N_r+1)*N_theta,N+1);
Fn = zeros((N_r+1)*N_theta,N+1);


%% default
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

[D,x] = cheb(N_r);

r = (d*x+b+a)/2;

% D_r = d/2*D;
% Dr = diag(r)*D_r;

%% Required 
A = diag(kron(f,b-r));
B = diag(kron(f_theta,b-r));
F = kron(diag(f),eye(N_r+1)); 
F_theta = kron(diag(f_theta),eye(N_r+1));
R = kron(eye(N_theta),diag(r));

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
xi = Ar*besselh(pp,k.*(a+Eps.*f)).*exp(1i*pp.*theta);
u_exact = Ar * kron(exp(1i*pp*theta),besselh(pp,k*r)); 
u_exact_theta = 1i*pp*Ar * kron(exp(1i*pp*theta),besselh(pp,k*r)); 
u_exact_r = Ar*k* kron(exp(1i*pp*theta),diff_besselh(pp,1,k*r));

%% tests
% u_r = op_Partial_r(u_exact,D,N_theta,d);
% norm(u_exact_r-u_r);
% 
% u_theta = op_D_theta(u_exact,N_r,N_theta,p);
% norm(u_exact_theta-u_theta);
Un_C=k^2*op_C(Un,A,F,R,a,b,N);
Fn(:,2)=(d*A*op_Partial_r(op_Dr(Un(:,1),D,N_theta,d,R),D,N_theta,d)+d*op_Dr(A*op_Partial_r(Un(:,1),D,N_theta,d),D,N_theta,d,R)...
    -d*F*op_D_theta(op_D_theta(Un(:,1),N_r,N_theta,p),N_r,N_theta,p)-d*op_D_theta(F*op_D_theta(Un(:,1),N_r,N_theta,p),N_r,N_theta,p)...
    -d*B*op_Partial_r(op_D_theta(Un(:,1),N_r,N_theta,p),D,N_theta,d)-d*op_D_theta(B*op_Partial_r(Un(:,1),D,N_theta,d),N_r,N_theta,p)...
    +d*F_theta*op_D_theta(Un(:,1),N_r,N_theta,p)+Un_C(:,1))/(-d^2);
for i=3:N+1
Fn(:,i)=(d*A*op_Partial_r(op_Dr(Un(:,i-1),D,N_theta,d,R),D,N_theta,d)+d*op_Dr(A*op_Partial_r(Un(:,i-1),D,N_theta,d),D,N_theta,d,R)...
    -d*F*op_D_theta(op_D_theta(Un(:,i-1),N_r,N_theta,p),N_r,N_theta,p)-d*op_D_theta(F*op_D_theta(Un(:,i-1),N_r,N_theta,p),N_r,N_theta,p)...
    -d*B*op_Partial_r(op_D_theta(Un(:,i-1),N_r,N_theta,p),D,N_theta,d)-d*op_D_theta(B*op_Partial_r(Un(:,i-1),D,N_theta,d),N_r,N_theta,p)...
    +d*F_theta*op_D_theta(Un(:,i-1),N_r,N_theta,p)...
    +A*op_Partial_r(A*op_Partial_r(Un(:,i-2),D,N_theta,d),D,N_theta,d)+F*op_D_theta(F*op_D_theta(Un(:,i-2),N_r,N_theta,p),N_r,N_theta,p)...
    +B*op_Partial_r(F*op_D_theta(Un(:,i-2),N_r,N_theta,p),D,N_theta,d)+F*op_D_theta(B*op_Partial_r(Un(:,i-2),D,N_theta,d),N_r,N_theta,p)...
    +B*op_Partial_r(B*op_Partial_r(Un(:,i-2),D,N_theta,d),D,N_theta,d)-F*F_theta*op_D_theta(Un(:,i-2),N_r,N_theta,p)...
    -F_theta*b*op_Partial_r(Un(:,i-2),D,N_theta,d)...
    +Un_C(:,i-2))/(-d^2);
end