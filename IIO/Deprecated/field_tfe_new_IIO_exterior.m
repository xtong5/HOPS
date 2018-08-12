function [Un,Dr_Un,Dp_Un] = field_tfe_new_IIO_exterior(xi_n,f,f_theta,k,a,b,p,N_theta,N,N_r,sigma,Y)

%% Default
d = b-a;
T_p = -k* diff_besselh(p,1,k*b)./besselh(p,k*b);

%% pre allocation
u_n = zeros(N_theta*(N_r+1),N+1);
u_n_p = zeros(N_theta*(N_r+1),N+1);
Dr_u_n = zeros((N_r+1)*N_theta,N+1);
Dp_u_n = zeros((N_r+1)*N_theta,N+1);
Un = zeros(N_theta,N+1);
Dr_Un = zeros(N_theta,N+1);
Dp_Un = zeros(N_theta,N+1);

%% data
[D,x] = cheb(N_r);
r = (d*x+b+a)/2;
Dr = D*(2/d); % partial_r 
A = diag(kron(f,b-r));
B = diag(kron(f_theta,b-r));
F = kron(diag(f),eye(N_r+1)); 
F_theta = kron(diag(f_theta),eye(N_r+1));
R = kron(eye(N_theta),diag(r));
C1 = -2*d*F*R*R + 2*d*A*R;
C2 = F*F*R*R - 4*F*A*R + A*A;
C3 = (2/d)*A*F*F*R - (2/d)*F*A*A;
C4 = (1/d^2)*F*F*A*A;

%% main loop
D_p = diff_mat_new_exterior(p,T_p,Dr,r,k,sigma,a,Y);
Zero = zeros(N_theta*(N_r+1),1);

for n=1:N+1
    
    if n == 1 % order=0
    xi_n_p = fft(xi_n(:,n));
    FF = zeros(N_r+1,1); 
    for i=1:N_theta
        rhs = FF;
        rhs(1) = 0;rhs(end) =xi_n_p(i);
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
    
    if n == 2 % order=1
    xi_n_p = fft(xi_n(:,n));
    K_n_p = K_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),Zero,Zero,Zero,Dr_u_n(:,n-1),Zero,...
        Dp_u_n(:,n-1),Zero,A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n); 
    Fn_p = FT(Fn,N_r,N_theta);
    Ln_p=L_n_new_exterior(u_n(:,n-1),Zero,xi_n(:,n-1),zeros(N_theta,1),Dr_u_n(:,n-1),...
        Zero,Dp_u_n(:,n-1),Zero,N_r,N_theta,f,f_theta,a,d,sigma,Y,n);
    for i=1:N_theta
        rhs = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        rhs(1) = K_n_p(i);rhs(end) =xi_n_p(i)+Ln_p(i);       
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n == 3 % order=2
    xi_n_p = fft(xi_n(:,n));
    K_n_p = K_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),u_n(:,n-2),Zero,Zero,Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n); 
    Fn_p = FT(Fn,N_r,N_theta);
    Ln_p=L_n_new_exterior(u_n(:,n-1),u_n(:,n-2),xi_n(:,n-1),xi_n(:,n-2),Dr_u_n(:,n-1),...
        Dr_u_n(:,n-2),Dp_u_n(:,n-1),Dp_u_n(:,n-2),N_r,N_theta,f,f_theta,a,d,sigma,Y,n);
    for i=1:N_theta
        rhs = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        rhs(1) = K_n_p(i);rhs(end) =xi_n_p(i)+Ln_p(i);
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n == 4 % order=3
    xi_n_p = fft(xi_n(:,n));
    K_n_p = K_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),u_n(:,n-2),u_n(:,n-3),Zero,Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
        Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    Ln_p=L_n_new_exterior(u_n(:,n-1),u_n(:,n-2),xi_n(:,n-1),xi_n(:,n-2),Dr_u_n(:,n-1),...
        Dr_u_n(:,n-2),Dp_u_n(:,n-1),Dp_u_n(:,n-2),N_r,N_theta,f,f_theta,a,d,sigma,Y,n);
    for i=1:N_theta
        rhs = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        rhs(1) = K_n_p(i);rhs(end) =xi_n_p(i)+Ln_p(i);
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end

    if n>4 % order=4 to N
    xi_n_p = fft(xi_n(:,n));
    K_n_p = K_n_exterior(u_n(:,n-1),N_r,N_theta,f,T_p,d);
    Fn=F_n_exterior(u_n(:,n-1),u_n(:,n-2),u_n(:,n-3),u_n(:,n-4),Dr_u_n(:,n-1),Dr_u_n(:,n-2),...
    Dp_u_n(:,n-1),Dp_u_n(:,n-2),A,B,D,F,F_theta,R,C1,C2,C3,C4,p,d,k,N_r,N_theta,n);
    Fn_p = FT(Fn,N_r,N_theta);
    Ln_p=L_n_new_exterior(u_n(:,n-1),u_n(:,n-2),xi_n(:,n-1),xi_n(:,n-2),Dr_u_n(:,n-1),...
        Dr_u_n(:,n-2),Dp_u_n(:,n-1),Dp_u_n(:,n-2),N_r,N_theta,f,f_theta,a,d,sigma,Y,n);
    for i=1:N_theta
        rhs = Fn_p((i-1)*(N_r+1)+1:i*(N_r+1));
        rhs(1) = K_n_p(i);rhs(end) =xi_n_p(i)+Ln_p(i);
        u_n_p((i-1)*(N_r+1)+1:i*(N_r+1),n)= D_p(:,:,i)\rhs;
    end
    u_n(:,n) = FT_inv(u_n_p(:,n),N_r,N_theta);
    Dr_u_n(:,n) = op_Partial_r(u_n(:,n),D,N_theta,d);
    Dp_u_n(:,n) = op_D_theta(u_n(:,n),N_r,N_theta,p);
    end
end


%% output
for i = 1:N_theta
    Un(i,:) = u_n(i*(N_r+1),:);
end

for i = 1:N_theta
    Dr_Un(i,:) = Dr_u_n(i*(N_r+1),:);
    Dp_Un(i,:) = Dp_u_n(i*(N_r+1),:);
end

end
