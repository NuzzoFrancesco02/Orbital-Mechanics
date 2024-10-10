function [t, E, M] = kepler_inv_equation(e, a, mu, t0, th)
% This function computes the time corresponding to a given true anomaly (th),
% taking into account multiple orbital periods.
%
% PROTOTYPE:
% [t, E, M] = kepler_inv_equation(e, a, mu, t0, th)
%
% INPUT:
% e [1x1]    - Orbital eccentricity [-]
% a [1x1]    - Semi-major axis (km) [L]
% mu [1x1]   - Gravitational parameter (km^3/s^2) [L^3/T^2]
% t0 [1x1]   - Initial time (s) [T]
% th [nx1]   - True anomaly (rad) [-]
%
% OUTPUT:
% t [nx1]    - Time corresponding to the true anomaly (s) [T]
% E [nx1]    - Eccentric anomaly (rad) [-]
% M [nx1]    - Mean anomaly (rad) [-]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

    % Compute eccentric anomaly (E) from the true anomaly (th)
    E = 2 * atan2(sqrt((1 - e) / (1 + e)) .* tan(th ./ 2), 1);

    % Mean motion (n) is computed as sqrt(mu / a^3)
    n = sqrt(mu / a^3);  % rad/s
    
    % Orbital period (T) is calculated as T = 2*pi / n
    T = 2 * pi / n;  % s

    % Ensure that the eccentric anomaly E is within the range [0, 2*pi]
    E = mod(E, 2 * pi);

    % Mean anomaly (M) is calculated from E using Kepler's equation
    M = E - e .* sin(E);
    
    % Compute the time corresponding to the mean anomaly in a single period
    t_single = t0 + M / n;  % Time in a single period
    
    % Determine the number of complete orbital periods by checking if th exceeds 2*pi
    num_periods = floor(th / (2 * pi));  % Number of completed orbital periods

    % Adjust the time by accounting for the number of complete periods
    t = t_single + num_periods * T;

end
