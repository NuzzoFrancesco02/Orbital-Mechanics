function [r, v, t, theta] = kep2car_time(t, a, e, i, OM, w, th0, mu)
% This function converts Keplerian orbital elements to position (r) and velocity (v) 
% vectors in the Earth-Centered Inertial (ECI) frame. It takes into account either 
% a single time instant or a time range. If the time vectore is composed by
% initial and final time, it's computationally more efficient
%
% PROTOTYPE:
% [r, v, t] = kep2ECI(t, a, e, i, OM, w, th0, mu)
%
% INPUT:
% t [nx1]    - Time vector (if two elements, calculates over a range) [T]
% a [1x1]    - Semi-major axis (km) [L]
% e [1x1]    - Eccentricity [-]
% i [1x1]    - Inclination (rad) [-]
% OM [1x1]   - Longitude of the ascending node (rad) [-]
% w [1x1]    - Argument of perigee (rad) [-]
% th0 [1x1]  - Initial true anomaly (rad) [-]
% mu [1x1]   - Gravitational parameter (km^3/s^2) [L^3/T^2]
%
% OUTPUT:
% r [nx3]    - Position vector in the ECI frame (rx, ry, rz) [L]
% v [nx3]    - Velocity vector in the ECI frame (vx, vy, vz) [L/T]
% t [nx1]    - Time vector corresponding to the computed positions (s) [T]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

    % Check if the time vector t has more than two elements
    if length(t) > 2
        % If there are more than two time values, calculate the true anomaly (theta)
        theta = elliptic.kepler_equation(t, e, a, mu, t(1), th0);  % Compute the true anomaly
        % Convert Keplerian elements to Cartesian coordinates (position and velocity)
        [r, v] = coord.kep2car_theta(a, e, i, OM, w, theta, mu);  % Convert to ECI

    elseif length(t) == 2
        % If there are exactly two time points, calculate for the range
        % Compute the true anomaly at the beginning and the end of the time range
        theta_beg = elliptic.kepler_equation(t(1), e, a, mu, t(1), th0);  % At t(1)
        theta_end = elliptic.kepler_equation(t(2), e, a, mu, t(1), th0);  % At t(2)
        
        % Interpolate the true anomaly values between the start and end points
        theta = linspace(theta_beg, theta_end, 1e4)';  % Linspace for smooth transition
        
        % Convert Keplerian elements to Cartesian coordinates over the interpolated range
        [r, v] = coord.kep2car_theta(a, e, i, OM, w, theta, mu);  % Convert to ECI
        
        % Recalculate the time vector using inverse Kepler's equation
        t = elliptic.kepler_inv_equation(e, a, mu, t(1), th0, theta);  % Time corresponding to true anomaly
    end

end
