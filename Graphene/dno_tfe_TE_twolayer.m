function [U_n,W_n] = dno_tfe_TE_twolayer(zeta_n,psi_n,f,f_theta,rho,...
    p,k_u,k_w,a,b,c,N_theta,N,N_r)

%% set up
U_n = zeros(N_theta,N+1);
W_n = zeros(N_theta,N+1);
Gn_u_U = zeros(N_theta,N+1,N+1); %dno for G_u_jU_i
Gn_w_W = zeros(N_theta,N+1,N+1); %dno for G_w_jW_i
Nn_W = zeros(N_theta,N+1); %normal operator for N_jW_i
G0_u = -a*k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
G0_w = a*k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);
Delta = 1./(G0_u + G0_w - rho);


%% compute order n=0
n = 0;
s = n;
R1hat = fft( zeta_n(:,0+1) );
R2hat = fft( -psi_n(:,0+1) );
U_n(:,0+1) = ifft((R1hat.*(G0_w-rho) + R2hat).*Delta);
W_n(:,0+1) = U_n(:,0+1) - zeta_n(:,0+1);

% compute and store Gn_u[U_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = U_n(:,s+1);
[~,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
Gn = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
for r=0:N-s
  Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
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
  R1hat = fft(zeta_n(:,n+1));
  R2hat = fft(-psi_n(:,n+1));
  for m=0:n-1
    R2hat = R2hat-fft(Gn_u_U(:,n+1-m,m+1))-fft(Gn_w_W(:,n+1-m,m+1))+rho.*fft(Nn_W(:,n+1-m,m+1));
  end
  U_n(:,n+1) = ifft((R1hat.*(G0_w-rho) + R2hat).*Delta);
  W_n(:,n+1) = U_n(:,n+1) - zeta_n(:,n+1);

  % compute and store Gn_u[U_s] 
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  [~,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
  Gn = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
  for r=0:N-s
    Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
  end
  
  % compute and store Gn_w[W_s] and |N|n[W_n]
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = W_n(:,s+1);
  [~,Dr_Wn,Dp_Wn] = field_tfe_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
  Gn = dno_tfe_helmholtz_polar_interior(Dr_Wn,Dp_Wn,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
  Nn = current_n(f_theta,xi,N-s);
  for r=0:N-s
      Gn_w_W(:,r+1,s+1) = Gn(:,r+1);
      Nn_W(:,r+1,s+1) = Nn(:,r+1);
  end

end