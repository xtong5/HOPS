function [Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_interior(xi_n,f,f_theta,k,a,c,p,N_theta,N,N_r)
%% Default
d=a-c;
alpha = (a+c)/d;
w = k*d/2;
T_p = k* diff_besselj(p,1,k*c)./besselj(p,k*c);


%% pre allocation
u_n = zeros(N_theta*(N_r+1),N+1);
u_n_p = zeros(N_theta*(N_r+1),N+1);
Dr_u_n = zeros((N_r+1)*N_theta,N+1);
Dp_u_n = zeros((N_r+1)*N_theta,N+1);
J_n_p = zeros(N_theta,N+1);
Un = zeros(N_theta,N+1);

%% data
[D,x] = cheb(N_r);
r = (d*x+c+a)/2;
A = diag(kron(f,r-c));
B = diag(kron(f_theta,r-c));
F = kron(diag(f),eye(N_r+1)); 
F_theta = kron(diag(f_theta),eye(N_r+1));
R = kron(eye(N_theta),diag(r));
X = (x+alpha).^2.*w^2;
C1 = 2*d*F*R*R + 2*d*A*R;
C2 = F*F*R*R + 4*F*A*R + A*A;
C3 = (2/d)*A*F*F*R + (2/d)*F*A*A;
C4 = (1/d^2)*F*F*A*A;


%% main loop
T_p_d = T_p*d/2;
D_p = diff_mat_interior(p,T_p_d,alpha,w,D,x);
Zero = zeros(N_theta*(N_r+1),1);
% xi_n_p = fft(xi_n);

for n=1:N+1    
    if n == 1 % order=0
    xi_n_p = fft(xi_n(:,n));
    J_n_p(:,n) = J_n_interior(u_n(:,n),N_r,N_theta,f,T_p,d);
    FF = zeros(N_r+1,1); 
    for i=1:N_theta
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i)).*x...
            +(-J_n_p(i,n)+xi_n_p(i)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i)).*(x+alpha)/(1+2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
    
    if n == 2 % order=1
    xi_n_p = fft(xi_n(:,n));
    J_n_p(:,n) = J_n_interior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_interior(u_n(:,n-1),Zero,Zero,Zero,Dr_u_n(:,n-1),Zero,...
        Dp_u_n(:,n-1),Zero,A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i)).*x...
            +(-J_n_p(i,n)+xi_n_p(i)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i)).*(x+alpha)/(1+2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
    
    if n == 3 % order=2
    xi_n_p = fft(xi_n(:,n));
    J_n_p(:,n) = J_n_interior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_interior(u_n(:,n-1),u_n(:,n-2),Zero,Zero,Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i)).*x...
            +(-J_n_p(i,n)+xi_n_p(i)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i)).*(x+alpha)/(1+2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n == 4 % order=3
    xi_n_p = fft(xi_n(:,n));
    J_n_p(:,n) = J_n_interior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_interior(u_n(:,n-1),u_n(:,n-2),u_n(:,n-3),Zero,Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i)).*x...
            +(-J_n_p(i,n)+xi_n_p(i)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i)).*(x+alpha)/(1+2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
    
    if n>4 % order=4 to N
    xi_n_p = fft(xi_n(:,n));
    J_n_p(:,n) = J_n_interior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_interior(u_n(:,n-1),u_n(:,n-2),u_n(:,n-3),u_n(:,n-4),Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    for i=1:N_theta
        FF = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        h_p = (J_n_p(i,n)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i)).*x...
            +(-J_n_p(i,n)+xi_n_p(i)+T_p_d(i)*xi_n_p(i))/(1+2*T_p_d(i));
        rhs = FF - (X-p(i)^2).*h_p-(J_n_p(i,n)+T_p_d(i)*xi_n_p(i)).*(x+alpha)/(1+2*T_p_d(i));
        rhs(1) = 0;rhs(end) =0;
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs+h_p;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
end

%% output
for i = 1:N_theta
    Un(i,:) = u_n(1+(i-1)*(N_r+1),:);
end

for i = 1:N_theta
    Dr_Un(i,:) = Dr_u_n(1+(i-1)*(N_r+1),:);
    Dp_Un(i,:) = Dp_u_n(1+(i-1)*(N_r+1),:);
end

end