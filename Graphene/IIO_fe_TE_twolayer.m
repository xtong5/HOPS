function [I_u_n,I_w_n] = IIO_fe_TE_twolayer(zeta_n,psi_n,f,f_theta,rho,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p)

%% set up
I_u_n = zeros(N_theta,N+1);
I_w_n = zeros(N_theta,N+1);
Nn_Iw = zeros(N_theta,N+1,N+1); %normal operator for N_jI_w_i
Qn = zeros(N_theta,N+1,N+1); %IIO for Q_u_jI_u_i
Sn = zeros(N_theta,N+1,N+1); %IIO for S_w_jI_w_i
Sn_sum= zeros(N_theta,N+1); %sum of S_w_jI_w_i
Normal_n = Nnorm_n(f_theta,N);

G0_u = -a*k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
G0_w = a*k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);
Q0_p = (G0_u + Z_p)./(G0_u + Y_p); S0_p = (G0_w - Y_p)./(G0_w - Z_p); 
Delta1 = S0_p - rho./(Y_p-Z_p).*(1-S0_p); Delta2 = 1 - rho./(Y_p-Z_p).*(1-S0_p);
Delta = 1./(Delta2 - Q0_p.*Delta1);

%% compute order n=0
n = 0;
s = n;
R1hat = Y_p.*fft(zeta_n(:,0+1)) - fft(psi_n(:,0+1));
R2hat = Z_p.*fft(zeta_n(:,0+1)) - fft(psi_n(:,0+1));
I_u_n(:,0+1) = ifft((Delta2.*R1hat - Delta1.*R2hat).*Delta);
I_w_n(:,0+1) = ifft((R2hat - Q0_p.*R1hat).*Delta);

% compute and store Qn[I_u_0]
xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_u_n(:,s+1); 
anp= field_fe_IIO_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Y_p);
Qn_fe = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Z_p);
for r=0:N-s
  Qn(:,r+1,s+1) = Qn_fe(:,r+1);
end

% compute and store Sn[I_w_0] and |N|n[I_w_0] and |N|nSn[I_w_0]
xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_w_n(:,s+1);
dnp = field_fe_IIO_helmholtz_polar_interior(xi,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Z_p);
Sn_fe = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Y_p);
Nn = current_n(f_theta,xi,N-s);
for r=0:N-s
  Sn(:,r+1,s+1) = Sn_fe(:,r+1);
  Nn_Iw(:,r+1,s+1) = Nn(:,r+1);
end

% compute and store sum{Sn[I_w_0]} 
Sn_sum(:,0+1) = Sn(:,0+1,0+1);


%% order n>0
for n=1:N
  s=n;
  R1hat = Y_p.*fft(zeta_n(:,n+1)) - fft(psi_n(:,n+1));
  R2hat = Z_p.*fft(zeta_n(:,n+1)) - fft(psi_n(:,n+1));
  for m=0:n-1
    R1hat = R1hat + fft( -(1+rho./(Y_p-Z_p)).*Sn(:,n-m+1,m+1) + rho./(Y_p-Z_p).*Nn_Iw(:,n-m+1,m+1))...
            - rho./(Y_p-Z_p).*fft( Normal_n(:,n-m+1).*Sn_sum(:,m+1));
    R2hat = R2hat + fft( -Qn(:,n-m+1,m+1) + rho./(Y_p-Z_p).*Nn_Iw(:,n-m+1,m+1) -rho./(Y_p-Z_p).*Sn(:,n-m+1,m+1))...
            - rho./(Y_p-Z_p).*fft( Normal_n(:,n-m+1).*Sn_sum(:,m+1));
  end
  I_u_n(:,n+1) = ifft((Delta2.*R1hat - Delta1.*R2hat).*Delta);
  I_w_n(:,n+1) = ifft((R2hat - Q0_p.*R1hat).*Delta);

  % compute and store Qn[I_u_s]
  xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_u_n(:,s+1); 
  anp= field_fe_IIO_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Y_p);
  Qn_fe = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Z_p);
  for r=0:N-s
      Qn(:,r+1,s+1) = Qn_fe(:,r+1);
  end
  
  % compute and store Sn[I_w_s] and |N|n[I_w_s] and |N|nSn[I_w_s]
  xi = zeros(N_theta,N-s+1); xi(:,0+1) = I_w_n(:,s+1);
  dnp = field_fe_IIO_helmholtz_polar_interior(xi,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Z_p);
  Sn_fe = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Y_p);
  Nn = current_n(f_theta,xi,N-s);
  for r=0:N-s
    Sn(:,r+1,s+1) = Sn_fe(:,r+1);
    Nn_Iw(:,r+1,s+1) = Nn(:,r+1);
  end
  
  % compute and store sum{Sn[I_w_l]} 
  for r=0:s
    Sn_sum(:,n+1) = Sn_sum(:,n+1)+Sn(:,r+1,s-r+1);
  end
  

end