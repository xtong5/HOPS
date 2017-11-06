function [D_p] = diff_mat_exterior(p,T_p,c,w,D,x)

%% pre allocation
n = size(x,1);
N_theta = size(p,1);
D_p = zeros(n,n,N_theta);
Diag = diag(x+c);
A = eye(n); 
DD = Diag*D*Diag*D + Diag*Diag*w^2;

for i=1:N_theta
    D_p(:,:,i) = DD - p(i)^2*A;
    D_p(end,:,i) = A(end,:);
    D_p(1,:,i)=D(1,:);
    D_p(1,1,i) = D_p(1,1,i)-T_p(i);
%     fprintf('Det = %d Cond = %d \n', det(D_p(:,:,i)),cond(D_p(:,:,i)));
end


end