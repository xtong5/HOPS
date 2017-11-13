% test_TFE_exterior.m
%
% Script to test Helmholtz DNO solvers in polar (exterior)
%
% XT 10/17


%% Default
k=2.1;
Eps = 0.02;
N_theta = 64;
N_r = 16;
a = 1;
b = 1.6;
d = b-a;
N = 16; 
L = 2*pi;
p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';
c = (b+a)/d;
w = k*d/2;
T_p = k* diff_besselh(p,1,k*b)./besselh(p,k*b);

fprintf('test_helmholtz_polar_TFE\n');
fprintf('-------------\n');
fprintf('k = %g a = %g b = %g\n',k,a,b);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');

%% pre allocation
u_n = zeros(N_theta*(N_r+1),N+1);
u_n_p = zeros(N_theta*(N_r+1),N+1);
Dr_u_n = zeros((N_r+1)*N_theta,N+1);
Dp_u_n = zeros((N_r+1)*N_theta,N+1);
J_n_p = zeros(N_theta,N+1);
Un = zeros(N_theta,N+1);
Dr_Un = zeros(N_theta,N+1);
Dp_Un = zeros(N_theta,N+1);
Gn_tfe = zeros(N_theta,N+1);

%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );
[D,x] = cheb(N_r);
r = (d*x+b+a)/2;
A = diag(kron(f,b-r));
B = diag(kron(f_theta,b-r));
F = kron(diag(f),eye(N_r+1)); 
F_theta = kron(diag(f_theta),eye(N_r+1));
R = kron(eye(N_theta),diag(r));
X = (x+c).^2.*w^2;
C1 = -2*d*F*R*R + 2*d*A*R;
C2 = F*F*R*R - 4*F*A*R + A*A;
C3 = (2/d)*A*F*F*R - (2/d)*F*A*A;
C4 = (1/d^2)*F*F*A*A;

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
xi = Ar*besselh(pp,k.*(a+Eps.*f)).*exp(1i*pp.*theta);
% u_exact = Ar * kron(exp(1i*pp*theta),besselh(pp,k*r));
u_exact = Ar * exp(1i*pp*theta) .* besselh(pp,k.*(a+Eps.*f));
% u_exact = Ar * exp(1i*pp*theta) .* besselh(pp,k.*(b));
xi_n_p = xi_n_exterior(N_theta,N,f,Ar,pp,a,k,theta);
nu =Ar*(-k.*(a+Eps.*f).*(diff_besselh(pp,1,k.*(a+Eps.*f)))+1i*pp*Eps.*...
    f_theta.*besselh(pp,k.*(a+Eps.*f))./(a+Eps.*f) ).*exp(1i*pp.*theta); %DNO

%% main loop
T_p_d = T_p*d/2;
D_p = diff_mat_exterior(p,T_p_d,c,w,D,x);
Zero = zeros(N_theta*(N_r+1),1);

for n=1:N+1
    
    if n == 1 % order=0
    J_n_p(:,n) = J_n_exterior(u_n(:,n),N_r,N_theta,f,T_p,d);
    FF = zeros(N_r+1,1); 
    for i=1:N_theta
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i)).*x...
            +(J_n_p(i,n)+xi_n_p(i,n)-T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n)).*(x+c)/(1-2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
    
    if n == 2 % order=1
    J_n_p(:,n) = J_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),Zero,Zero,Zero,Dr_u_n(:,n-1),Zero,...
        Dp_u_n(:,n-1),Zero,A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n); 
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i)).*x...
            +(J_n_p(i,n)+xi_n_p(i,n)-T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n)).*(x+c)/(1-2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n == 3 % order=2
    J_n_p(:,n) = J_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),u_n(:,n-2),Zero,Zero,Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n); 
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i)).*x...
            +(J_n_p(i,n)+xi_n_p(i,n)-T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n)).*(x+c)/(1-2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n == 4 % order=3
    J_n_p(:,n) = J_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),u_n(:,n-2),u_n(:,n-3),Zero,Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i)).*x...
            +(J_n_p(i,n)+xi_n_p(i,n)-T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n)).*(x+c)/(1-2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n>4 % order=4 to N
    J_n_p(:,n) = J_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),u_n(:,n-2),u_n(:,n-3),u_n(:,n-4),Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
    Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i)).*x...
            +(J_n_p(i,n)+xi_n_p(i,n)-T_p_d(i)*xi_n_p(i,n))/(1-2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i,n)).*(x+c)/(1-2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
end


%% compare 
for i = 1:N_theta
    Un(i,:) = u_n(i*(N_r+1),:);
%     Un_p(i,:) = u_n_p(i*(N_r+1),:);
%     Un(i,:) = u_n(1+(i-1)*(N_r+1),:);
end
% [relerr,nplot] = compute_errors_2d_polar(u_exact,Un,Eps,N,N_theta);

for i = 1:N_theta
    Dr_Un(i,:) = Dr_u_n(i*(N_r+1),:);
    Dp_Un(i,:) = Dp_u_n(i*(N_r+1),:);
end
Gn_tfe(:,1) = -a*Dr_Un(:,1);
Gn_tfe(:,2) = -f.*(1/a-1/d).*Gn_tfe(:,1)-a*Dr_Un(:,2)-2*f.*Dr_Un(:,1)...
    +f_theta.*Dp_Un(:,1)/a;
for i=3:N+1
    Gn_tfe(:,i) = -f.*(1/a-1/d).*Gn_tfe(:,i-1)+f.*f.*Gn_tfe(:,i-2)/(a*d)...
        -a*Dr_Un(:,i)-2*f.*Dr_Un(:,i-1)-(f.*f+f_theta.*f_theta).*Dr_Un(:,i-2)/a...
        +f_theta.*Dp_Un(:,i-1)/a-f_theta.*f.*Dp_Un(:,i-2)/(a*d);
end
[relerr,nplot] = compute_errors_2d_polar(nu,Gn_tfe,Eps,N,N_theta);        

