function [th, E, M] = kepler_equation(t, e, a, mu, t0, th0, toll)
% kepler_equation - Solves Kepler's equation for the true anomaly, eccentric anomaly, and mean anomaly.
%
% PROTOTYPE:
% [th, E, M] = kepler_equation(t, e, a, mu, t0, th0, toll)
%
% INPUT:
% t [nx1]    - Time values from the initial time [T]
% e [1]      - Eccentricity of the orbit [-]
% a [1]      - Semi-major axis of the orbit [L]
% mu [1]     - Gravitational parameter [L^3/T^2]
% t0 [1]     - Initial time [T]
% th0 [1]    - Initial true anomaly at time t0 [rad]
% toll [1]   - Tolerance for the numerical solver [default: 1e-6]
%
% OUTPUT:
% th [nx1]   - True anomaly [rad]
% E [nx1]    - Eccentric anomaly [rad]
% M [nx1]    - Mean anomaly [rad]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

    % Set default tolerance if not provided
    if nargin == 6
        toll = 1e-6;
    end
    
    % Calculate mean motion
    n = sqrt(mu / a^3); T = 2*pi/n;
    n_orbs = floor(t./T);
    % Initial eccentric anomaly (E0)
    E0 = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(th0 / 2));
    
    % Initial guess for eccentric anomaly (E)
    dt = t - t0;
    E_guess = n * dt + (e * sin(n * dt)) ./ (1 - sin(n * dt + e) + sin(n * dt));

    % Set options for the numerical solver
    options = optimoptions('fsolve', 'FunctionTolerance', toll, 'Display', 'none');
    
    % Solve Kepler's equation for eccentric anomaly (E)
    E = fsolve(@(E) n * dt - E + e * sin(E) + E0 - e * sin(E0), E_guess, options);
    
    % Calculate true anomaly (th) from eccentric anomaly (E)
    th = 2 * atan2(sqrt((1 + e) / (1 - e)) * tan(E / 2), 1) + n_orbs.*2*pi;
    
    % Calculate mean anomaly (M)
    M = n * dt;
end
