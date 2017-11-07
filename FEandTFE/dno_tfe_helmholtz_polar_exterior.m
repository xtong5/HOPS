function [Gn_tfe] = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k,a,b,p,N_theta,N,N_r)
d=b-a;
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

end