
  
m=1;
for j=1:N_theta
    k = floor(N/2);
    Gu(j,m) = padesum(Gn_fe_u(j,:).',Eps,k);    
end
U(m) = norm(Gu(:,m),2)/sqrt(N_theta);      

