Mode = 2; 
L = 2*pi;
k=2.1;

sigma = 0;
% eta = 1.5;
eta = 1;

% Eps = 0.02;
Eps = 0;
N_theta = 64;
N = 16;
N_r = 16;
a = 1.6;
c = 0.8;
d = a-c;

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';
T_p = k * diff_besselj(p,1,k*c)./besselj(p,k*c);

[D,x] = cheb(N_r);
r = (d*x+c+a)/2;
Dr = D*(2/d); % partial_r 


D_p1 = zeros(N_r+1,N_r+1,N_theta);
D_p = zeros(N_r+1,N_r+1,N_theta);
Diag = diag(r);
DD = Diag*Dr*Diag*Dr + k*k*Diag*Diag;

for i=1:N_theta
    D_p1(:,:,i) = DD - p(i)^2*A;
    D_p1(1,:,i) = A(1,:);
    D_p1(end,:,i)=Dr(end,:);
    D_p1(end,end,i) = D_p1(end,end,i)-T_p(i);
%     fprintf('Det = %d Cond = %d \n', det(D_p(:,:,i)),cond(D_p(:,:,i)));
end


for i=1:N_theta
    D_p(:,:,i) = DD - p(i)^2*A;
%     D_p(1,:,i) = 1/sigma*a*Dr(1,:);
    D_p(1,:,i) = sigma*a*Dr(1,:);
%     D_p(1,1,i) = D_p(1,1,i) + 1i*eta;
    D_p(1,1,i) = D_p(1,1,i) + eta;
    D_p(end,:,i) = Dr(end,:);
    D_p(end,end,i) = D_p(end,end,i)-T_p(i);
%     fprintf('Det = %d Cond = %d \n', det(D_p(:,:,i)),cond(D_p(:,:,i)));
end

norm(D_p1(:,:,1)-D_p(:,:,1))

