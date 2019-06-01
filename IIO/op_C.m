function [Fn_C]=op_C(u_nmo,u_nmt,u_nmth,u_nmf,C1,C2,C3,C4,n)
% terms by C1,C2,C3,C4
if n == 1
    Fn_C = zeros(size(u_nmo));
end
if n == 2
    Fn_C = C1*u_nmo;
end
if n == 3
    Fn_C = C1*u_nmo+C2*u_nmt;
end
if n == 4
    Fn_C = C1*u_nmo+C2*u_nmt+C3*u_nmth; 
end
if n > 4
    Fn_C = C1*u_nmo+C2*u_nmt+C3*u_nmth+C4*u_nmf;
end
end