function [anp] = field_fe_helmholtz_polar_exterior1(xi,f,k,p,N_theta,N,a)
anp = zeros(N_theta,N+1); 
fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
anp(:,0+1) = fft(xi); 
for n=1:N
  anp(:,n+1) = field_fe_n_helmholtz_polar_exterior(xi,fn,k,p,n,a,anp);  
end