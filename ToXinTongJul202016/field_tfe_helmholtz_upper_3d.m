function [un] = field_tfe_helmholtz_upper_3d(xi_n,f,p1,p2,alpha1p,alpha2p,betap,...
    eep,eem,Dy,a,Nx1,Nx2,Ny,N)

un = zeros(Nx1,Nx2,Ny+1,N+1);
xihat_n = zeros(Nx1,Nx2,N+1);

k2 = alpha1p(0+1)^2 + alpha2p(0+1)^2 + betap(0+1)^2;

ell_top = 0 + 1;
for n=0:N
  xihat_n(:,:,n+1) = fft2(eem.*xi_n(:,:,n+1));
end
f_x1 = real(ifft2( (1i*p1).*fft2(f) ));
f_x2 = real(ifft2( (1i*p2).*fft2(f) ));

ll = [0:Ny]';
y_min = 0.0; y_max = a;
y = ((y_max-y_min)/2.0)*(cos(pi*ll/Ny) - 1.0) + y_max;

f_full = zeros(Nx1,Nx2,Ny+1);
f_x1_full = zeros(Nx1,Nx2,Ny+1);
f_x2_full = zeros(Nx1,Nx2,Ny+1);
for ell=0:Ny
  f_full(:,:,ell+1) = f;
  f_x1_full(:,:,ell+1) = f_x1;
  f_x2_full(:,:,ell+1) = f_x2;
end
Uhat = zeros(Nx1,Nx2,Ny+1);

a_minus_y_full = zeros(Nx1,Nx2,Ny+1);
for j1=1:Nx1
  for j2=1:Nx2
    for ell=0:Ny
      a_minus_y_full(j1,j2,ell+1) = a - y(ell+1);
    end
  end
end

% Order zero

for ell=0:Ny
  un(:,:,ell+1,0+1) = eep.*ifft2( exp(1i*betap*y(ell+1)).*xihat_n(:,:,0+1) );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx1,Nx2,Ny+1);
  Jn = zeros(Nx1,Nx2);
  
  A1_x1x1 = -(2.0/a)*f_full;
  %A1_x1x2 = 0
  A1_x1y = -(1.0/a)*(a_minus_y_full).*f_x1_full;
  %A1_x2x1 = 0
  A1_x2x2 = -(2.0/a)*f_full;
  A1_x2y = -(1.0/a)*(a_minus_y_full).*f_x2_full;
  A1_yx1 = A1_x1y;
  A1_yx2 = A1_x2y;
  %A1_yy = 0
  
  A2_x1x1 = (1.0/a^2)*f_full.^2;
  %A2_x1x2 = 0
  A2_x1y = (1.0/a^2)*(a_minus_y_full).*(f_full.*f_x1_full);
  %A2_x2x1 = 0
  A2_x2x2 = (1.0/a^2)*f_full.^2;
  A2_x2y = (1.0/a^2)*(a_minus_y_full).*(f_full.*f_x2_full);
  A2_yx1 = A2_x1y;
  A2_yx2 = A2_x2y;
  A2_yy = (1.0/a^2)*((a_minus_y_full).^2).*(f_x1_full.^2 + f_x2_full.^2);
  
  B1_x1 = -(1.0/a)*f_x1_full;
  B1_x2 = -(1.0/a)*f_x2_full;
  %B1_y = 0
  
  B2_x1 = (1.0/a^2)*f_full.*f_x1_full;
  B2_x2 = (1.0/a^2)*f_full.*f_x2_full;
  B2_y = (1.0/a^2).*(a_minus_y_full).*(f_x1_full.^2 + f_x2_full.^2);
  
  C1 = -k2*(2.0/a)*f_full;
  C2 = k2*(1.0/a^2)*f_full.^2;
  
  if(n>=1)
    u_x1 = dx1p_3d(un(:,:,:,n-1+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_x1x1.*u_x1;
    Fn = Fn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    %A1_x2x1 = 0
    temp = A1_yx1.*u_x1;
    Fn = Fn - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B1_x1.*u_x1;
    Fn = Fn + temp;
    
    u_x2 = dx2p_3d(un(:,:,:,n-1+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A1_x1x2 = 0
    temp = A1_x2x2.*u_x2;
    Fn = Fn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_yx2.*u_x2;
    Fn = Fn - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B1_x2.*u_x2;
    Fn = Fn + temp;
    
    u_y = dy_3d(un(:,:,:,n-1+1),Dy,a,Nx1,Nx2,Ny);
    temp = A1_x1y.*u_y;
    Fn = Fn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_x2y.*u_y;
    Fn = Fn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A1_yy = 0
    %B1_y = 0
    
    temp = C1.*un(:,:,:,n-1+1);
    Fn = Fn - temp;
    
    Su = eep.*ifft2( (1i*betap).*fft2(eem.*un(:,:,ell_top,n-1+1)) );
    Jn = Jn - (1.0/a)*f.*Su;
  end
  
  if(n>=2)
    u_x1 = dx1p_3d(un(:,:,:,n-2+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x1x1.*u_x1;
    Fn = Fn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    %A2_x2x1 = 0
    temp = A2_yx1.*u_x1;
    Fn = Fn - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B2_x1.*u_x1;
    Fn = Fn + temp;
    
    %A2_x1x2 = 0
    u_x2 = dx2p_3d(un(:,:,:,n-2+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x2x2.*u_x2;
    Fn = Fn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_yx2.*u_x2;
    Fn = Fn - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B2_x2.*u_x2;
    Fn = Fn + temp;
    
    u_y = dy_3d(un(:,:,:,n-2+1),Dy,a,Nx1,Nx2,Ny);
    temp = A2_x1y.*u_y;
    Fn = Fn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x2y.*u_y;
    Fn = Fn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_yy.*u_y;
    Fn = Fn - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B2_y.*u_y;
    Fn = Fn + temp;
    
    temp = C2.*un(:,:,:,n-2+1);
    Fn = Fn - temp;
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx1,Nx2,Ny+1);
  for ell=0:Ny
    Fnhat(:,:,ell+1) = fft2(eem.*Fn(:,:,ell+1));
  end
  Jnhat = fft2(eem.*Jn);
  Fnhat_p = zeros(Ny+1,1);
  
  for j1=1:Nx1
    for j2=1:Nx2
      for ell=0:Ny
        Fnhat_p(ell+1) = Fnhat(j1,j2,ell+1);
      end
      alphaalpha = 1.0;
      betabeta = 0.0;
      gammagamma = k2 - alpha1p(j1,j2)^2 - alpha2p(j1,j2)^2;
      d_a = 1.0;
      n_a = 0.0;
      r_a = xihat_n(j1,j2,n+1);
      d_b = -1i*betap(j1,j2);
      n_b = 1.0;
      r_b = Jnhat(j1,j2);
      %{
      uhat_p = solvebvp(Fnhat_p,alphaalpha,betabeta,gammagamma,y_min,y_max,...
          d_a,n_a,r_a,d_b,n_b,r_b);
      %}
      %{
      uhat_p = solvebvp_alt(Fnhat_p,alphaalpha,betabeta,gammagamma,y_min,y_max,...
          d_a,n_a,r_a,d_b,n_b,r_b);
      %}
      uhat_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
          (2.0/(y_max-y_min))*Dy,d_a,n_a,r_a,d_b,n_b,r_b);
      
      for ell=0:Ny
        Uhat(j1,j2,ell+1) = uhat_p(ell+1);
      end
    end
  end
  
  for ell=0:Ny
    un(:,:,ell+1,n+1) = eep.*ifft2(Uhat(:,:,ell+1));
  end

end