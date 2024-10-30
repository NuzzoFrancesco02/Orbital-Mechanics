function [r, v, theta] = kep2car_theta(a, e, i, OM, om, theta, mu)
% kep2car_theta - Conversion from Keplerian elements to Cartesian coordinates given theta.
%                 Handles both elliptical (e < 1), parabolic (e = 1), and hyperbolic (e > 1) orbits.
%                 Supports vectorized input for the true anomaly [1xN].
%
% PROTOTYPE:
%   [r, v, theta] = kep2car_theta(a, e, i, OM, om, theta, mu)
%
% INPUT:
%   a      [1x1]   - Semi-major axis                         [km]
%   e      [1x1]   - Eccentricity                            [-]
%   i      [1x1]   - Inclination                             [rad]
%   OM     [1x1]   - Right Ascension of the Ascending Node   [rad]
%   om     [1x1]   - Argument of Pericenter                  [rad]
%   theta  [1xN]   - True anomaly                            [rad]
%   mu     [1x1]   - Gravitational parameter                 [km^3/s^2]
%
% OUTPUT:
%   r      [3xN]   - Position vector in Cartesian coordinates [km]
%   v      [3xN]   - Velocity vector in Cartesian coordinates [km/s]
%   theta  [1xN]   - True anomaly (useful in cases where e >= 1) [rad]
%
% CONTRIBUTORS:
%   Francesco Nuzzo
%
% VERSIONS:
%   2024-10-10: First version
%
% -------------------------------------------------------------------------

% Calculate the semi-latus rectum of the orbit
p = a * (1 - e^2);

% If the orbit is parabolic or hyperbolic, adjust theta range
if e == 1  % Parabolic case
    p = 2 * a;  
    rp_max = 1e5;  % Arbitrary maximum radius to limit theta range
    if length(theta) ~= 1
        theta = -floor(acos(p / rp_max - 1)):0.005:floor(acos(p / rp_max - 1));
    end
elseif e > 1  % Hyperbolic case
    rp_max = 1e5;  
    if length(theta) ~= 1
        theta = -floor(acos((p / rp_max - 1) / e)):0.01:floor(acos((p / rp_max - 1) / e));
    end
end


% Transformation matrix from perifocal to Cartesian coordinates
R3 = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R1 = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
R = R1 * R2 * R3;
R = R';

% Initialize position and velocity vectors
r = []; 
v = [];

% Loop through each angle in theta
for j = 1:length(theta)
    % Calculate position and velocity in perifocal coordinates
    if e < 1  % Elliptical orbit
        r_pf = (p / (1 + e * cos(theta(j)))) * [cos(theta(j)); sin(theta(j)); 0];
        v_pf = sqrt(mu / p) * [-sin(theta(j)); e + cos(theta(j)); 0];
    elseif e > 1  % Hyperbolic orbit
        r_pf = real((p / (1 + e * cos(theta(j)))) * [cos(theta(j)); sin(theta(j)); 0]);
        v_pf = imag(sqrt(mu / p) * [-sin(theta(j)); e + cos(theta(j)); 0]);
    else  % Parabolic orbit (e == 1)
        r_pf = (p / (1 + e * cos(theta(j)))) * [cos(theta(j)); sin(theta(j)); 0];
        v_pf = sqrt(mu / p) * [-sin(theta(j)); e + cos(theta(j)); 0];
    end
    
    % Convert position and velocity to Cartesian coordinates
    r = [r, R * r_pf];
    
    v = [v, R * v_pf];
  
end
