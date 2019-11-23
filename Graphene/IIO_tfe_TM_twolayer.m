function [I_u_n,I_w_n] = IIO_tfe_TM_twolayer(Nzeta_n,psi_n,f,f_theta,eta,...
    p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,Y_p,Z_p)

%% set up
I_u_n = zeros(N_theta,N+1);
I_w_n = zeros(N_theta,N+1);
Nn_Iu = zeros(N_theta,N+1,N+1); %normal operator for N_jI_u_i
Nn_Iw = zeros(N_theta,N+1,N+1); %normal operator for N_jI_w_i
Qn = zeros(N_theta,N+1,N+1); %IIO for Q_u_jI_u_i
Sn = zeros(N_theta,N+1,N+1); %IIO for S_w_jI_w_i
Qn_sum= zeros(N_theta,N+1); %sum of Q_u_jI_u_i
Sn_sum= zeros(N_theta,N+1); %sum of S_w_jI_w_i

Npsi_n = current_n(f_theta,psi_n,N);
Normal_n = Nnorm_n(f_theta,N);

A = -sigma_u*a*k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
B = sigma_w*a*k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);
Q0_p = (A + Z_p)./(A + Y_p);
S0_p = (B - Y_p)./(B - Z_p); 
Delta1 = S0_p + eta.*Y_p./(2.*sigma_w).*(1+S0_p); 
Delta2 = 1 + eta.*Z_p./(2.*sigma_w).*(1+S0_p);
Delta = 1./(Delta2 - Q0_p.*Delta1);

%% compute order n=0
n = 0;
s = n;
R1hat = Y_p.*fft(Nzeta_n(:,0+1)) - sigma_u.*fft(Npsi_n(:,0+1));
R2hat = Z_p.*fft(Nzeta_n(:,0+1)) - sigma_u.*fft(Npsi_n(:,0+1));
I_u_n(:,0+1) = ifft((Delta2.*R1hat - Delta1.*R2hat).*Delta );
I_w_n(:,0+1) = ifft( (R2hat - Q0_p.*R1hat).*Delta );

% compute and store Qn[I_u_0] and |N|n[I_u_0]
xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_u_n(:,s+1); 
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r,sigma_u,Y_p);
Qn_tfe = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N-s,sigma_u,Z_p);
Nn = current_n(f_theta,xi,N-s);
for r=0:N-s
  Qn(:,r+1,s+1) = Qn_tfe(:,r+1);
  Nn_Iu(:,r+1,s+1) = Nn(:,r+1);
end

% compute and store Sn[I_w_0] and |N|n[I_w_0]
xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_w_n(:,s+1);
[Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r,sigma_w,Z_p);
Sn_tfe = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N-s,sigma_w,Y_p);
Nn = current_n(f_theta,xi,N-s);
for r=0:N-s
  Sn(:,r+1,s+1) = Sn_tfe(:,r+1);
  Nn_Iw(:,r+1,s+1) = Nn(:,r+1);
end

% compute and store sum{Sn[I_u_0]} and sum{Sn[I_w_0]} 
Qn_sum(:,0+1) = Qn(:,0+1,0+1);
Sn_sum(:,0+1) = Sn(:,0+1,0+1);


%% order n>0
for n=1:N
  s=n;
  R1hat = Y_p.*fft(Nzeta_n(:,n+1)) - sigma_u.*fft(Npsi_n(:,n+1));
  R2hat = Z_p.*fft(Nzeta_n(:,n+1)) - sigma_u.*fft(Npsi_n(:,n+1));
  for m=0:n-1
    R1hat = R1hat + fft( -Nn_Iu(:,n-m+1,m+1) - Sn(:,n-m+1,m+1) ...
        - 0.5*eta.*Y_p./sigma_w.*Sn(:,n-m+1,m+1) - Normal_n(:,n-m+1).*Sn_sum(:,m+1));
    R2hat = R2hat + fft( -Qn(:,n-m+1,m+1) - Normal_n(:,n-m+1).*Qn_sum(:,m+1) ...
        - Nn_Iw(:,n-m+1,m+1) - 0.5*eta.*Z_p./sigma_w.*Sn(:,n-m+1,m+1));
  end
  I_u_n(:,n+1) = ifft((Delta2.*R1hat - Delta1.*R2hat).*Delta);
  I_w_n(:,n+1) = ifft((R2hat - Q0_p.*R1hat).*Delta);

  % compute and store Qn[I_u_s] and |N|n[I_u_s] and |N|nSn[I_u_s]
  xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_u_n(:,s+1); 
  [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r,sigma_u,Y_p);
  Qn_tfe = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N-s,sigma_u,Z_p);
  Nn = current_n(f_theta,xi,N-s);
  for r=0:N-s
      Qn(:,r+1,s+1) = Qn_tfe(:,r+1);
      Nn_Iu(:,r+1,s+1) = Nn(:,r+1);
  end
  
  % compute and store Sn[I_w_s] and |N|n[I_w_s] and |N|nSn[I_w_s]
  xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_w_n(:,s+1);
  [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r,sigma_w,Z_p);
  Sn_tfe = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N-s,sigma_w,Y_p);
  Nn = current_n(f_theta,xi,N-s);
  for r=0:N-s
    Sn(:,r+1,s+1) = Sn_tfe(:,r+1);
    Nn_Iw(:,r+1,s+1) = Nn(:,r+1);
  end
  
  % compute and store sum{Qn[I_u_l]} and sum{Sn[I_w_l]} 
  for r=0:s
    Qn_sum(:,n+1) = Qn_sum(:,n+1)+Qn(:,r+1,s-r+1);
    Sn_sum(:,n+1) = Sn_sum(:,n+1)+Sn(:,r+1,s-r+1);
  end
  

end