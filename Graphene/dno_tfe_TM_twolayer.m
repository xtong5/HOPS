function [U_n,W_n] = dno_tfe_TM_twolayer(Nzeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,eta,tau2,a,b,c,N_theta,N,N_r)

%% set up
U_n = zeros(N_theta,N+1);
W_n = zeros(N_theta,N+1);
Gn_u_U = zeros(N_theta,N+1,N+1); %dno for G_u_jU_i
Gn_w_W = zeros(N_theta,N+1,N+1); %dno for G_w_jU_i
Nn_U = zeros(N_theta,N+1); %normal operator for N_jU_i
Nn_W = zeros(N_theta,N+1); %normal operator for N_jW_i
% Nn_zeta = zeros(N_theta,N+1); %normal operator for N_jzeta_i

G0_u = -a*k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
G0_w = a*k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);
Delta = 1./(tau2.*G0_w+G0_u.*(1-eta.*G0_w));


% %% Compute and store |N|n[zeta]
% for n=0:N
%   s = n;
%   xi = zeros(N_theta,N-s+1);
%   xi(:,0+1) = zeta_n(:,s+1);
%   Nn = current_n(f_theta,xi,N-s);
%   for r=0:N-s
%     Nn_zeta(:,r+1,s+1) = Nn(:,r+1);
%   end
% end


%% compute order n=0
n = 0;
s = n;
R1hat = fft(Nzeta_n(:,0+1));
R2hat = fft(-psi_n(:,0+1));
U_n(:,0+1) = ifft((R1hat.*tau2.*G0_w + R2hat.*(1-eta.*G0_w)).*Delta);
W_n(:,0+1) = ifft((R2hat-R1hat.*G0_u).*Delta);

% compute and store Gn_u[U_0] and |N|n[U_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = U_n(:,s+1);
[~,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
Gn = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
Nn = current_n(f_theta,xi,N-s);
for r=0:N-s
  Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
  Nn_U(:,r+1,s+1) = Nn(:,r+1);
end

% compute and store Gn_w[W_0] and |N|n[W_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = W_n(:,s+1);
[~,Dr_Wn,Dp_Wn] = field_tfe_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
Gn = dno_tfe_helmholtz_polar_interior(Dr_Wn,Dp_Wn,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
Nn = current_n(f_theta,xi,N-s);
for r=0:N-s
  Gn_w_W(:,r+1,s+1) = Gn(:,r+1);
  Nn_W(:,r+1,s+1) = Nn(:,r+1);
end


%% order n>0
for n=1:N
  s=n;
%   R1hat = 0;
  R1hat = fft(Nzeta_n(:,n+1));
  R2hat = -fft(psi_n(:,n+1));
  for m=0:n-1
    R1hat = R1hat + fft(-Nn_U(:,n-m+1,m+1)+Nn_W(:,n-m+1,m+1)-eta.*Gn_w_W(:,n-m+1,m+1));
    R2hat = R2hat - fft(Gn_u_U(:,n-m+1,m+1)) - tau2.*fft(Gn_w_W(:,n-m+1,m+1));
  end
%   for m=0:n
%     R1hat = R1hat + fft(Nn_zeta(:,n-m+1,m+1));    
%   end
  U_n(:,n+1) = ifft((R1hat.*tau2.*G0_w + R2hat.*(1-eta.*G0_w)).*Delta);
  W_n(:,n+1) = ifft((R2hat-R1hat.*G0_u).*Delta);
  
  % compute and store Gn_u[U_s] and |N|n[U_s]
  xi = zeros(N_theta,N-s+1); xi(:,0+1) = U_n(:,s+1); 
  [~,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
  Gn = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
  Nn = current_n(f_theta,xi,N-s);
  for r=0:N-s
    Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
    Nn_U(:,r+1,s+1) = Nn(:,r+1);
  end
  
  % compute and store Gn_w[W_s] and |N|n[W_s]
  xi = zeros(N_theta,N-s+1); xi(:,0+1) = W_n(:,s+1);
  [~,Dr_Wn,Dp_Wn] = field_tfe_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
  Gn = dno_tfe_helmholtz_polar_interior(Dr_Wn,Dp_Wn,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
  Nn = current_n(f_theta,xi,N-s);
  for r=0:N-s
    Gn_w_W(:,r+1,s+1) = Gn(:,r+1);
    Nn_W(:,r+1,s+1) = Nn(:,r+1);
  end

end