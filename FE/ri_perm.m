function [n,epsilon] = ri_perm(lambda,material)

n = 0.0;
epsilon = 0.0;

lambda_nm = 1000*lambda;

if(strcmp(material,'VACUUM'))
  n = 1.0;
  epsilon = n^2;
elseif(strcmp(material,'EPOXY'))
  n = 1.5375...
      + 8290.45/(lambda_nm^2.0)...
      - 2.11046e8/(lambda_nm^4.0);
  epsilon = n^2;
elseif(strcmp(material,'GOLD'))
  % Lorentz Model: Values for Gold
  % Comment: Lambda_Min=0.27, Lambda_Max=6
  % Rakic, Djurisic, Elazar, & Majewski, Appl. Opt. 37, 5271 (1998)
  
  omega = 2.0*pi/lambda;

  epsilon_infty = 1.0;

  Delta = [1589.516;50.19525;20.91469;148.4943;1256.973;9169.0];
  a = [1.0;1.0;1.0;1.0;1.0;1.0];
  b = [0.268419;1.220548;1.747258;4.406129;12.63;11.21284];
  c = [0.0;4.417455;17.66982;226.0978;475.1387;4550.765];
  
  epsilon = epsilon_infty;

  for j=1:6
    epsilon = epsilon + Delta(j)/(-a(j)*omega^2.0-1i*b(j)*omega+c(j));
  end

  n = sqrt(epsilon);
elseif(strcmp(material,'WATER'))
  n2 = 1.0...
      + (5.666959820e-1)*(lambda^2.0)/((lambda^2.0)-5.084151894e-3)...
      + (1.731900098e-1)*(lambda^2.0)/((lambda^2.0)-1.818488474e-2)...
      + (2.095951857e-2)*(lambda^2.0)/((lambda^2.0)-2.625439472e-2)...
      + (1.125228406e-1)*(lambda^2.0)/((lambda^2.0)-1.073842352e1);
  n = sqrt(n2);
  epsilon = n2;
elseif(strcmp(material,'SILVER'))
  omega = 2.0*pi/lambda;
  
  epsilon_infty = 1.0;

  Delta = [1759.471;135.344;258.1946;22.904036;1749.06;11756.18];
  a = [1.0;1.0;1.0;1.0;1.0;1.0];
  b = [0.243097;19.68071;2.289161;0.329194;4.639097;12.25];
  c = [0.0;17.07876;515.022;1718.3357;2116.092;10559.42];

  epsilon = epsilon_infty;

  for j=1:6
    epsilon = epsilon + Delta(j)/(-a(j)*(omega^2.0)-1i*b(j)*omega+c(j));
  end
  
  n = sqrt(epsilon);
elseif(strcmp(material,'TUNGSTEN'))
  % Tungsten (Lorentz Model)
  % 3/23/15: DPN Interpretation
  %
  % Lorentz Model: Values for Tungsten
  % Comment: Lambda_Min=0.081, Lambda_Max=4.03
  %          0.1 eV - 5 eV
  % Source: App.Optics/Vol-37/pp-5271
  
  % Note: E (eV) = 1.23984187/lambda (microns)
  %       omega = 2*pi/lambda
  %       E = (1.23984187/(2*pi)) omega = (1/sigma) omega
  
  sigma = (2*pi)/1.23984187;

  omega = 2.0*pi/lambda;
  
  epsilon_infty = 1.0;

  omegap = 13.22;
    
  fj = [0.206;0.054;0.166;0.706;2.590];
  Gammaj = [0.064;0.530;1.281;3.332;5.836];
  omegaj = [0;1.004;1.917;3.580;7.498];
  
  for j=1:5
    Delta(j) = fj(j)*(omegap*sigma)^2.0;
    a(j) = 1.0;
    b(j) = Gammaj(j)*sigma;
    c(j) = (omegaj(j)*sigma)^2.0;
  end
  
  epsilon = epsilon_infty;

  for j=1:5
    epsilon = epsilon + Delta(j)/(-a(j)*(omega^2.0)-1i*b(j)*omega+c(j));
  end
  
  n = sqrt(epsilon);
end

return;