function [anp] = field_fe_helmholtz(xi,f,p,alphap,gammap,eep,eem,Nx,N)
anp = zeros(Nx,N+1); fn = zeros(Nx,N+1);
fn(:,0+1) = ones(Nx,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
anp(:,0+1) = fft( eem.*xi ); 
for n=1:N
  anp(:,n+1) = 0*xi;
  for m=0:n-1
    anp(:,n+1) = anp(:,n+1) - fft(eem.*( fn(:,n-m+1).*...
        eep.*ifft((1i*gammap).^(n-m).*anp(:,m+1)) )); %??
  end
end