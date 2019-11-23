function [Un_far] = B_fe_farfield(anp,k,a,b,p,N_theta,N)
Un_far = zeros(N_theta,N+1); 

for n=1:N+1
    Un_far(:,n)=ifft(anp(:,n).*besselh(p,k.*b)./besselh(p,k.*a)); 
end


end