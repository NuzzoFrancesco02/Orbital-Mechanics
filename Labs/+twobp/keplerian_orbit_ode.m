function [t, r, v] = keplerian_orbit_ode(y0, tspan, mu, str, J2, Re)
% keplerian_orbit_ode Propagates the orbit of a celestial body under 
% the influence of gravitational forces using either a two-body or J2 model.
%
% PROTOTYPE
% [t, r, v] = keplerian_orbit_ode(y0, tspan, mu, str, J2, Re)
%
% INPUT:
% y0[6x1]       Initial state of the body (rx, ry, rz, vx, vy, vz) [ L, L/T ]
% tspan[nx1]    Time interval for the simulation [T]
% mu[1]         Gravitational parameter of the primary body [L^3/T^2]
% str           Optional string; set to 'J2' to include the J2 perturbation model.
% J2[1]         Gravitational field constant associated with the J2 effect [-]
% Re[1]         Mean radius of the planet [L]
%
% OUTPUT:
% t[nx1]       Time vector for each time step [T]
% r[nx1]       Position vector for each time step [L]
% v[nx1]       Velocity vector for each time step [L/T]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS
% 2024-10-08: First version
%
% -------------------------------------------------------------------------

% Initial conditions: Extract position and velocity from the input vector.
r0_norm = norm(y0(1:3)); % Magnitude of the initial position vector.
v0_norm = norm(y0(4:6)); % Magnitude of the initial velocity vector.

% Compute angular momentum vector and eccentricity vector.
h_0 = cross(y0(1:3), y0(4:6)); % Specific angular momentum.
e_0 = cross(y0(4:6), h_0) / mu - y0(1:3) / r0_norm; % Eccentricity vector.

% Plot the energy profile for the initial conditions.
plots.energy_plot(y0(1:3), y0(4:6), mu);

%% ODE Integration
% Set options for the ODE solver with high precision.
options = odeset('RelTol', 1e-15, 'AbsTol', 1e-15);

% Choose the appropriate ODE function based on the presence of J2 model.
if nargin < 4
    % Call the two-body problem ODE function if J2 is not specified.
    [t, y] = ode113(@(t, y) twobp.ode_2bp(t, y, mu), tspan, y0, options);
elseif nargin == 6 && strcmp(str,'J2')
    % Call the two-body problem with J2 perturbation if specified.
    [t, y] = ode113(@(t, y) twobp.ode_2bp_J2(t, y, mu, J2, Re), tspan, y0, options);
else
    % Display an error message for insufficient input parameters.
    fprintf('\nAttention! Insufficient number of inputs');
end

% Separate the solution into position and velocity vectors.
r = y(:, 1:3); % Position vector at each time step.
v = y(:, 4:6); % Velocity vector at each time step.

%% Solution Check
% Analyze the accuracy of the computed solution.
plots.sol_errors(t, r, v, mu);

%% Velocity Components
% Plot velocity components as a function of time.
plots.velocity_components(t, r, v);

%% Plot the Orbit
% Visualize the orbit in 3D space.
plots.static_orbit(r);
