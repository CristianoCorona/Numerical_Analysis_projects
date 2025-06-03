% This script defines the main physical constants that are used in the DD mode
T      = 300;                              % ref. temperature value ( [] = K) 
kB     = 1.3806488e-23 ;                   % Boltzmann constant ([] = V A s K^{-1} = m^2 kg s^{-2} K^{-1})
q      = 1.6021892e-19 ;                   % elementary charge ([] = C)
Vth    = kB*T/q ;                          % thermal voltage ([] = V)
Nav    = 6.022e23 ;                        % Avogadro's constant ([] = mol^{-1})
eps0   = 8.854187818e-12 ;                 % vacuum dielectric permittivity ([] = F/m = A s V^{-1} m^{-1})
erw    = 80;                               % relative dielectric permittivity of water
erSi   = 11.7;                             % relative dielectric permittivity of silicon
erSiO2 = 3.9;                              % relative dielectric permittivity of silicon dioxide
m0     = 9.109534e-31;                     % electron rest mass ([] = Kg = V A s^3 m^{-2})
hPl    = 6.626176e-34;                     % Planck constant ([] = V A s^2 = J s)
hPlred = hPl/(2*pi);                       % reduced Planck constant ([] = V A s^2 = J s)
%
eps_Si = eps0*erSi;                        % dielectric constant of Silicon ([] F/m)
eps_Ox = eps0*erSiO2;                      % dielectric constant of Silicon Dioxide ([] F/m)
%
