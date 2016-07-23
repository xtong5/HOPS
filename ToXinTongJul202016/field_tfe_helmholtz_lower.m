function [un] = field_tfe_helmholtz_lower(xi_n,f,p,alphap,betap,eep,eem,Dy,b,Nx,Ny,N)

un = zeros(Nx,Ny+1,N+1);
xihat_n = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + betap(0+1)^2;

ell_bottom = Ny + 1;
for n=0:N
  xihat_n(:,n+1) = fft(eem.*xi_n(:,n+1));
end
f_x = real(ifft( (1i*p).*fft(f) ));

ll = [0:Ny]';
y_min = -b; y_max = 0.0;
y = ((y_max-y_min)/2.0)*(cos(pi*ll/Ny) - 1.0) + y_max;

f_full = zeros(Nx,Ny+1);
f_x_full = zeros(Nx,Ny+1);
for ell=0:Ny
  f_full(:,ell+1) = f;
  f_x_full(:,ell+1) = f_x;
end
Uhat = zeros(Nx,Ny+1);

y_plus_b_full = zeros(Nx,Ny+1);
for j=1:Nx
  y_plus_b_full(j,:) = y.' + b;
end

% Order zero

for ell=0:Ny
  un(:,ell+1,0+1) = eep.*ifft( exp(-1i*betap*y(ell+1)).*xihat_n(:,0+1) );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Ny+1);
  Jn = zeros(Nx,1);
  
  A1_xx = (2.0/b)*f_full;
  A1_xy = -(1.0/b)*(y_plus_b_full).*f_x_full;
  A1_yx = A1_xy;
  %A1_yy = 0;
  
  A2_xx = (1.0/b^2)*f_full.^2;
  A2_xy = -(1.0/b^2)*(y_plus_b_full).*(f_full.*f_x_full);
  A2_yx = A2_xy;
  A2_yy = (1.0/b^2)*((y_plus_b_full).^2).*(f_x_full.^2);
  
  B1_x = (1.0/b)*f_x_full;
  %B1_y = 0;
  
  B2_x = (1.0/b^2)*f_full.*f_x_full;
  B2_y = -(1.0/b^2).*(y_plus_b_full).*(f_x_full.^2);
  
  C1 = k2*(2.0/b)*f_full;
  C2 = k2*(1.0/b^2)*f_full.^2;
  
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Ny);
    temp = A1_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A1_yx.*u_x;
    Fn = Fn - dy(temp,Dy,b,Nx,Ny);
    temp = B1_x.*u_x;
    Fn = Fn + temp;
    
    u_y = dy(un(:,:,n-1+1),Dy,b,Nx,Ny);
    temp = A1_xy.*u_y;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Ny);
    %A1_yy = 0
    %B1_y = 0
    
    temp = C1.*un(:,:,n-1+1);
    Fn = Fn - temp;
    
    Su = eep.*ifft( (-1i*betap).*fft(eem.*un(:,ell_bottom,n-1+1)) );
    Jn = Jn + (1.0/b)*f.*Su;
  end
  
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Ny);
    temp = A2_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A2_yx.*u_x;
    Fn = Fn - dy(temp,Dy,b,Nx,Ny);
    temp = B2_x.*u_x;
    Fn = Fn + temp;
    
    u_y = dy(un(:,:,n-2+1),Dy,b,Nx,Ny);
    temp = A2_xy.*u_y;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A2_yy.*u_y;
    Fn = Fn - dy(temp,Dy,b,Nx,Ny);
    temp = B2_y.*u_y;
    Fn = Fn + temp;
    
    temp = C2.*un(:,:,n-2+1);
    Fn = Fn - temp;
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx,Ny+1);
  for ell=0:Ny
    Fnhat(:,ell+1) = fft(eem.*Fn(:,ell+1));
  end
  Jnhat = fft(eem.*Jn);
  
  for j=1:Nx
    Fnhat_p = Fnhat(j,:).';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - (alphap(j))^2;
    d_a = -(-1i*betap(j));
    n_a = 1.0;
    r_a = Jnhat(j);
    d_b = 1.0;
    n_b = 0.0;
    r_b = xihat_n(j,n+1);
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
    
    Uhat(j,:) = uhat_p.';
  end
  
  for ell=0:Ny
    un(:,ell+1,n+1) = eep.*ifft(Uhat(:,ell+1));
  end

end