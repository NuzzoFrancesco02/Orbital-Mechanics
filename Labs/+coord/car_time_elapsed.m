function [r, v, t] = car_time_elapsed(tspan, r0, v0, mu)
% This function computes the position and velocity in the Earth-Centered Inertial (ECI) 
% frame by numerically solving the two-body problem using the initial Cartesian 
% position and velocity vectors.
%
% PROTOTYPE:
% [r, v, t] = car_time_elapsed(tspan, r0, v0, mu)
%
% INPUT:
% tspan [nx1] - Time vector over which the integration is performed [T]
% r0 [3x1]    - Initial position vector (rx, ry, rz) [L]
% v0 [3x1]    - Initial velocity vector (vx, vy, vz) [L/T]
% mu [1x1]    - Gravitational parameter of the primary body (km^3/s^2) [L^3/T^2]
%
% OUTPUT:
% r [nx3]     - Position vectors over time in the ECI frame [L]
% v [nx3]     - Velocity vectors over time in the ECI frame [L/T]
% t [nx1]     - Time vector corresponding to the computed positions and velocities [T]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

    % Set the tolerances for the ODE solver (ode113) for higher accuracy
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

    % Ensure the initial position vector (r0) is a column vector
    if size(r0, 2) > 1
        r0 = r0';
    end
    
    % Ensure the initial velocity vector (v0) is a column vector
    if size(v0, 2) > 1
        v0 = v0';
    end
    
    % Combine the initial position and velocity into a single state vector
    s_0 = [r0; v0];  % State vector: position + velocity
    
    % Numerically integrate the two-body problem using ode113
    [t, y] = ode113(@(t, y) twobp.ode_2bp(t, y, mu), tspan, s_0, options);
    
    % Extract the position (r) and velocity (v) from the solution
    r = y(:, 1:3);  % Position vectors over time
    v = y(:, 4:6);  % Velocity vectors over time

end
