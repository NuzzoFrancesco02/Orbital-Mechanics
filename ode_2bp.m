function dy = ode_2bp(~, y, mu)
% ODE system for the two-body problem (Keplerian motion).
% 
% This function computes the derivatives of the state vector for a body
% under the influence of a central gravitational force (two-body problem).
%
% PROTOTYPE:
% dy = ode_2bp(t, y, mu)
%
% INPUT:
% t [1]        - Time (can be omitted, as the system is autonomous) [T]
% y [6x1]      - State vector: position (rx, ry, rz) and velocity (vx, vy, vz) [L, L/T]
% mu [1]       - Gravitational parameter of the primary body [L^3/T^2]
%
% OUTPUT:
% dy [6x1]     - Derivatives of the state: velocity and acceleration [L/T^2, L/T^3]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

    % Check if the input y is a column vector
    if size(y, 2) ~= 1
        error('\nAttention! y must be a column vector.');
    end

    % Extract the position and velocity vectors from the state vector y
    r = y(1:3);  % Position vector [rx, ry, rz] in [L]
    v = y(4:6);  % Velocity vector [vx, vy, vz] in [L/T]

    % Compute the norm (magnitude) of the position vector (distance from the primary body)
    r_norm = norm(r);  % Scalar distance [L]

    % The output dy contains:
    % 1. The velocity (derivative of the position) -> v
    % 2. The gravitational acceleration (Keplerian) -> -mu/r_norm^3 * r
    dy = [v;                   % Derivative of position is velocity
         (-mu / r_norm^3) * r]; % Gravitational acceleration in [L/T^2]
    
end
