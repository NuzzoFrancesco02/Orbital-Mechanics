function dy = ode_2bp_J2(t, y, mu, J2, Re)
% ODE system for the two-body problem (Keplerian motion) with J2 effect.
%
% This function computes the derivatives of the state vector for a satellite
% under the influence of Earth's gravitational field, including the J2 effect.
%
% PROTOTYPE:
% dy = ode_2bp_J2(t, y, mu, J2, Re)
%
% INPUT:
% t [1]         - Time (can be omitted, as the system is autonomous) [T]
% y [6x1]       - State vector: position (rx, ry, rz) and velocity (vx, vy, vz) [L, L/T]
% mu [1]        - Gravitational parameter of the primary body [L^3/T^2]
% J2 [1]        - J2 gravitational field constant [-]
% Re [1]        - Mean radius of the planet [L]
%
% OUTPUT:
% dy [6x1]      - Derivatives of the state: velocity and acceleration [L/T^2, L/T^3]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

    % Check that the input y is a column vector
    if size(y, 2) ~= 1
        error('\nAttention! y must be a column vector.');
    end

    % Extract position and velocity vectors from the state vector y
    r = y(1:3);  % Position vector [rx, ry, rz]
    v = y(4:6);  % Velocity vector [vx, vy, vz]
    
    % Compute the norm of the position vector
    r_norm = norm(r);  % Distance from the center of the primary body

    % Initialize the output dy with velocity and the gravitational acceleration without J2
    dy = [v;           % Derivative of position is velocity
         (-mu / r_norm^3) * r];  % Gravitational acceleration (Keplerian)

    % Calculate the J2 perturbation acceleration
    aj2 = (3/2) * J2 * mu * Re^2 / r_norm^4 * ...
          [0;  % No J2 contribution in x and y for a spherical primary
           0;
           0;
           r(1) / r_norm * (5 * r(3)^2 / r_norm^2 - 1);  % x-component of J2 acceleration
           r(2) / r_norm * (5 * r(3)^2 / r_norm^2 - 1);  % y-component of J2 acceleration
           r(3) / r_norm * (5 * r(3)^2 / r_norm^2 - 3)]; % z-component of J2 acceleration

    % Add the J2 perturbation to the derivative of the velocity (i.e., to the acceleration)
    dy = dy + aj2;

end
