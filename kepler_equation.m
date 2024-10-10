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
        toll = 1e-6; % Default tolerance for the solver
    end
    
    % Calculate elapsed time from the initial time
    t = t - t0; % Shift time to start from t0
    
    % Calculate mean motion and orbital period
    n = sqrt(mu / a^3); % Mean motion [rad/s]
    T = 2 * pi / n;     % Orbital period [s]
    
    % Determine the number of completed orbits
    k = floor(t ./ T); % Number of completed orbits
    
    % Initial guess for eccentric anomaly (E0)
    E0 = 2 * atan2(sqrt((1 - e) / (1 + e)) * tan(th0 ./ 2), 1);
    
    % Guess for the eccentric anomaly based on time
    E_guess = n .* t + (e .* sin(n .* t)) ./ (1 - sin(n .* t + e) + sin(n .* t));
    
    % Set options for the numerical solver
    options = optimoptions('fsolve', 'FunctionTolerance', toll);
    
    % Adjust time to the last completed orbit
    t = t - k .* T; % Adjusting time for the last complete orbit
    
    % Solve Kepler's equation for eccentric anomaly (E)
    E = fsolve(@(E) n .* t - E + e .* sin(E), E_guess, options) + E0;
    
    % Adjust E for completed orbits
    E = k * 2 * pi + E - E0; % Total eccentric anomaly adjusted for full orbits
    
    % Calculate true anomaly (th) from eccentric anomaly (E)
    th = 2 * atan2(sqrt((1 + e) / (1 - e)) * tan(E ./ 2), 1) + 2 * k * pi;
    
    % Calculate mean anomaly (M)
    M = n .* t; % Mean anomaly calculation
end
