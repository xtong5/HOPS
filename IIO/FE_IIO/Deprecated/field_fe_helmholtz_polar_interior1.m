function [dnp] = field_fe_helmholtz_polar_interior1(xi,f,k,p,N_theta,N,a)
dnp = zeros(N_theta,N+1); 
fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
dnp(:,0+1) = fft(xi); 
for n=1:N
  dnp(:,n+1) = field_fe_n_helmholtz_polar_interior(xi,fn,k,p,n,a,dnp);  
end
end