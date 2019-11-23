function [I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p)

%% set up
I_u_n = zeros(N_theta,N+1);
I_w_n = zeros(N_theta,N+1);
Qn = zeros(N_theta,N+1,N+1); %IIO for Q_u_jI_u_i
Sn = zeros(N_theta,N+1,N+1); %IIO for S_w_jI_w_i
A = -sigma_u*a*k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
B = sigma_w*a*k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);
Q0_p = (A + Z_p)./(A + Y_p);
S0_p = (B - Y_p)./(B - Z_p); 
Delta = 1./(1-Q0_p.*S0_p);



%% compute order n=0
n = 0;
s = n;
R1hat = Y_p.*fft(zeta_n(:,0+1))-sigma_u.*fft(psi_n(:,0+1));
R2hat = Z_p.*fft(zeta_n(:,0+1))-sigma_u.*fft(psi_n(:,0+1));
I_u_n(:,0+1) = ifft((R1hat - S0_p.*R2hat).*Delta);
I_w_n(:,0+1) = ifft((-Q0_p.*R1hat + R2hat).*Delta);

% compute and store Qn[I_u_0]
xi_u = zeros(N_theta,N-s+1); xi_u(:,0+1) = I_u_n(:,s+1); 
anp= field_IIO_fe_helmholtz_polar_exterior(xi_u,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Y_p);
Qn_fe = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Z_p);
for r=0:N-s
  Qn(:,r+1,s+1) = Qn_fe(:,r+1);
end

% compute and store Sn[I_w_0]
xi_w = zeros(N_theta,N-s+1); xi_w(:,0+1) = I_w_n(:,s+1);
dnp = field_IIO_fe_helmholtz_polar_interior(xi_w,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Z_p);
Sn_fe = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Y_p);
for r=0:N-s
  Sn(:,r+1,s+1) = Sn_fe(:,r+1);
end


%% order n>0
for n=1:N
  s=n;
  R1hat = Y_p.*fft(zeta_n(:,n+1))-sigma_u.*fft(psi_n(:,n+1));
  R2hat = Z_p.*fft(zeta_n(:,n+1))-sigma_u.*fft(psi_n(:,n+1));
  for m=0:n-1
    R1hat = R1hat - fft(Sn(:,n+1-m,m+1));
    R2hat = R2hat - fft(Qn(:,n+1-m,m+1));
  end
  I_u_n(:,n+1) = ifft((R1hat - S0_p.*R2hat).*Delta);
  I_w_n(:,n+1) = ifft((-Q0_p.*R1hat + R2hat).*Delta);

  % compute and store Qn[I_u_n]
  xi_u = zeros(N_theta,N-s+1); xi_u(:,0+1) = I_u_n(:,s+1); 
  anp= field_IIO_fe_helmholtz_polar_exterior(xi_u,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Y_p);
  Qn_fe = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s,sigma_u,Z_p);
  for r=0:N-s
    Qn(:,r+1,s+1) = Qn_fe(:,r+1);
  end
  
  % compute and store Sn[I_w_n]  
  xi_w = zeros(N_theta,N-s+1); xi_w(:,0+1) = I_w_n(:,s+1);
  dnp = field_IIO_fe_helmholtz_polar_interior(xi_w,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Z_p);
  Sn_fe = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s,sigma_w,Y_p);
  for r=0:N-s
    Sn(:,r+1,s+1) = Sn_fe(:,r+1);
  end

end