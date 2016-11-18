% this function computes and saves the data to file
function [lambda_crit,U_norm,W_norm,BU_norm] = ...
    FE_app_taylor(M,f,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name)

lambda = linspace(0.3,0.8,M);
epsvec = linspace(0,Eps_max,N_eps);
k_zero = (2*pi./lambda).';
n_u = zeros(M,1);
n_w = zeros(M,1);
epsilon_u = zeros(M,1);
epsilon_w = zeros(M,1);
epsilon_u_plot = zeros(M,1);
epsilon_w_plot = zeros(M,1);

for i = 1:M
    [n_u(i),epsilon_u(i)] = ri_perm(lambda(i),OUT);
    [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),IN);
    epsilon_u_plot(i) = epsilon_u(i);
    epsilon_w_plot(i) = epsilon_w(i);
end

[mini,r] = min(abs(epsilon_u - real(-epsilon_w)));
lambda_crit = lambda(r);

k_u = n_u.*k_zero; 
k_w = n_w.*k_zero; 
Mode = 2; 
% Mode = 1;


p = [0:N_theta/2-1,-N_theta/2:-1]';
f_theta = ifft(1i*p.*fft(f));

P_0 = 1;

BU = zeros(N_theta,M);
U = zeros(N_theta,M);
W = zeros(N_theta,M);
B_far = zeros(N_theta,N+1);
BU_norm = zeros(N_eps,M);
U_norm = zeros(N_eps,M);
W_norm = zeros(N_eps,M);
U_n = zeros(N_theta,N+1,M);
W_n = zeros(N_theta,N+1,M); 

for m =1:M
zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u(m)*a.*Sin);
zeta_n(:,0+1) = -Exp; 
psi_n(:,0+1) = (1i*k_u(m))*a*Sin.*Exp;

for n=1:N
    f_n = f.*f_n/n;
    if n > 1
        f_nmo = f.*f_nmo/(n-1);
    end
    zeta_n(:,n+1) = -Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
    psi_n(:,n+1) = (1i*k_u(m)).*Exp.*(a*(-1i*k_u(m))^n.*f_n.*...
    Sin.^(n+1)+(f.*Sin-f_theta.*cos(theta)).*...
    (-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
 end
if(Mode==1)
    tau2 = 1;
else
    tau2 = k_u(m)^2/k_w(m)^2;
end

   U_n(:,:,m) = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u(m),k_w(m),a,N_theta,N);
   apn_fe = field_fe_helmholtz_polar_exterior(U_n(:,:,m),f,k_u(m),a,p,N_theta,N);
   W_n(:,:,m) = U_n(:,:,m) - zeta_n;
   
   for n=1:N+1
     B_far(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u(m).*b)./besselh(p,k_u(m).*a));
   end 

   for l=1:N_eps
       Eps = epsvec(l);
       for j=1:N_theta
        BU(j,m) = taylorsum(B_far(j,:).',Eps,N);
        U(j,m) = taylorsum(U_n(j,:,m).',Eps,N);
        W(j,m) = taylorsum(W_n(j,:,m).',Eps,N);
       end
       
       U_norm(l,m) = norm(U(:,m),2)/sqrt(N_theta);
       BU_norm(l,m) = norm(BU(:,m),2)/sqrt(N_theta);
       W_norm(l,m) = norm(W(:,m),2)/sqrt(N_theta);
   end
   
end

filename = sprintf('data_%s.mat',name);
save(filename,'M','lambda','U_norm','W_norm','BU_norm',...
    'lambda_crit','epsilon_u_plot','epsilon_w_plot','epsvec','a','b',...
    'N_theta','N','N_eps','Eps_max','OUT','IN')




 

