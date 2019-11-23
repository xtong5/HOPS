function [Qn_tfe] = IIO_tfe_helmholtz_polar_exterior1(I_u_n,f,f_theta,k,a,b,p,N_theta,N,N_r,sigma,eta)
% should be I_u not I_u_n
d=b-a;
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k,a,b,p,N_theta,N,N_r,sigma,eta);
Gn_tfe = zeros(size(Un));

for i=1:N+1
    if i == 1
        Gn_tfe(:,1) = -a*Dr_Un(:,1);
    end
    if i == 2
        Gn_tfe(:,2) = -f.*(1/a-1/d).*Gn_tfe(:,1)-a*Dr_Un(:,2)-...
            2*f.*Dr_Un(:,1)+f_theta.*Dp_Un(:,1)/a;
    end
    if i > 2
        Gn_tfe(:,i) = -f.*(1/a-1/d).*Gn_tfe(:,i-1)+f.*f.*Gn_tfe(:,i-2)/(a*d)...
        -a*Dr_Un(:,i)-2*f.*Dr_Un(:,i-1)-(f.*f+f_theta.*f_theta).*Dr_Un(:,i-2)/a...
        +f_theta.*Dp_Un(:,i-1)/a-f_theta.*f.*Dp_Un(:,i-2)/(a*d);
    end
end

Qn_tfe = 1/sigma*Gn_tfe - 1i*eta*Un; 

end
