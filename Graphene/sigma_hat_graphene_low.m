function [sigma_hat_alma,sigma_hat_zsdlz,sigma_hat_drude,...
    sigma_hat_kubo,sigma_hat_drude_bfpv]...
    = sigma_hat_graphene_low(omega,T,mu,eta,E_F,Gamma)
% Dimensionless surface conductivity for graphene (in Siemens):
%
% \sigma = \epsilon_0 c_0 \sigma_hat
%
% Models from:
%
% a.) Angelis, Locatelli, Mutti, Aceves (ALMA) Opt. Lett. 41(3) 2016
%
% b.) Zhan, Shi, Dai, Liu, Zi (ZKDLZ) J. Phys. Condens. Matt. 25 2013
%
% c.) Drude (Tony Low, personal communication, August 2016)
%
% d.) Kubo (Tony Low, personal communication, August 2016)
%
% e.) Drude (Bludov, Ferreira, Peres, Vasilevskiy (BFPV) 
%     Int. J. Mod. Phys. B 27 (2013)
%
% Inputs:
%
% omega [rad/s]: frequency
% T [K]: temperature
% mu [eV]: chemical potential
% eta []: loss parameter (1 \leq eta \leq 5); DIMENSIONLESS

t = 2.7;                       % [eV]: hopping parameter
hbar = 6.582119514e-16;        % [eV s]: reduced Planck's constant
k_B = 8.6173324e-5;            % [eV/K]: Boltzmann constant
alpha = 0.0072973525664;       % []: Fine structure constant

Omega = hbar*omega/mu;            % []: Normalized frequency
kappa = 4.0*k_B*T/mu;          % []: Normalized temperature

sigma_hat_alma = 1i*(4.0*alpha/Omega)*(1.0 - (2.0/9.0)*(mu/t)^2)...
    + pi*alpha*eta*0.5*( tanh((Omega+2.0)/kappa) ...
    + tanh((Omega-2.0)/kappa) )...
    + 1i*alpha*log(abs((Omega-2.0)/(Omega+2.0)));

sigma_hat_zsdlz = 1i*(4.0*alpha/Omega)...
    + pi*alpha*heaviside_new(Omega-2.0)...
    + 1i*alpha*log(abs((Omega-2.0)/(Omega+2.0)));

% Drude and Kubo Models

%csp = 3e8;      % DPN: speed of light
hb = 1.055e-34; % DPN: hbar
q0 = 1.602e-19; % DPN: the electron charge? J/eV?
kb = 1.38e-23;  % DPN: Boltzman constant?
eta_low = 10e-3*q0; %hb/1e-12; % 1 ps
mu_low = mu*q0;
epsilon_0 = 8.854187817*1e-12;
c_0 = 299792458;
w = omega;

sigma_drude = 1i*q0^2*mu_low/(pi*hb^2*(w+1i*eta_low/hb));
sigma_hat_drude = sigma_drude/(epsilon_0*c_0);

zarray = linspace(0,18,100001)*q0/hb;
tol = 1e-3*q0/hb;
H1 = sinh(hb*zarray/(kb*T))./(cosh(mu_low/(kb*T))...
    +cosh(hb*zarray/(kb*T)));
H2 = sinh(hb*w/(2*kb*T))./(cosh(mu_low/(kb*T))+cosh(hb*w/(2*kb*T)));
Z = (H1-H2)./(w^2-4*zarray.^2);
Z(abs(2*zarray-w)<tol) = -hb/(kb*T)*cosh(mu_low/(kb*T))...
    *cosh(hb*w/2/(kb*T))/(cosh(mu_low/(kb*T))...
    +cosh(hb*w/2/(kb*T)))^2/(4*w);
sigma2 = q0^2*(H2+4*1i*w*trapz(zarray,Z)/pi)/4/hb;

den = pi*hb^2*(w+1i*eta_low/hb);
num = 1i*2*q0^2*kb*T*log(abs(2*cosh(mu_low/(2*kb*T))));
sigma_kubo = num/den+sigma2;
sigma_hat_kubo = sigma_kubo/(epsilon_0*c_0);

% Bludov, Ferreira, Peres, and Vasilevskiy Model (page 7, (23))

sigma_0 = pi*epsilon_0*c_0*alpha;
%E_F = 0.45;                       % 0.45 eV
%Gamma = 2.6*(1e-3);               % 2.6 meV
sigma_drude_bfpv = sigma_0*(4*E_F/pi)/(Gamma-1i*hbar*omega);
sigma_hat_drude_bfpv = sigma_drude_bfpv/(epsilon_0*c_0);

return;

function [hh] = heaviside_new(x)

if(x>0)
  hh = 1.0;
elseif(x==0)
  hh = 0.5;
else
  hh = 0.0;
end

return;