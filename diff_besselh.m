function [diffbesselh]=diff_besselh(p,n,z)
% nth derivative of hankel at 
% wavenumber vector p and scalar z OR
% wavebumber p and vector z
if size(z,1)==1
    diffbesselh=zeros(size(p,1),1);
    for i=1:n
        diffbesselh=diffbesselh+(-1)^i*nchoosek(n,i).*besselh(p'-n+2*i,z)';   
    end
end
if size(p,1)==1
    diffbesselh=zeros(size(z,1),1);
    for i=1:n
        diffbesselh=diffbesselh+(-1)^i*nchoosek(n,i).*besselh(p-n+2*i,z')';
    end
end
diffbesselh=diffbesselh./2^n;
end

