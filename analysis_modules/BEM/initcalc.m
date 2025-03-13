% perform initial calculations
function [Rp Sp q rps rpm r w TK TK0 mu0 mu] = initcalc(rR,Dp,rho,U0,J,TC)

Rp  = Dp/2;          % propeller radius          [m]
Sp  = (pi/4)*(Dp^2); % propeller disk area       [m2]
q   = .5*rho*(U0^2); % dynamic pressure          [Pa]
rps = U0/(J*Dp);     % revolutions per second    [1/s]
rpm = rps*60;        % revolutions per minute    [1/min]
r   = rR*Rp;         % dimensionfull radial coordinate [m]
w   = 2*pi*rps;      % rotational speed          [rad/s]
TK0 = 273.15;        % degrees Kelvin at 0 degrees Celsius
TK  = TC + TK0;      % air temperature           [K]
mu0 = 1.7894e-5;     % dynamic viscosity at sea level [kg/m*s]
mu  = mu0 * ((TK/TK0)^1.5) * ((TK0+110)/(TK+110)); % Sutherland formula for dynamic viscosity




