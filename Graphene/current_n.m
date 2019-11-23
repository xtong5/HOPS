function Js_n = current_n(f_theta,mu_n,N)
% norm of normal vector N operates on mu_n

N_theta = size(f_theta,1);
Js_n = zeros(N_theta,N+1);

Normal_n = Nnorm_n(f_theta,N);

for n=0:N   
    for m=0:n
        Js_n(:,n+1) = Js_n(:,n+1)+ Normal_n(:,n-m+1).*mu_n(:,m+1);
    end
end
