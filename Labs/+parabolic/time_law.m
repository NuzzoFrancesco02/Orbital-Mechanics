function [th, D, M] = time_law(t, rp, mu, t0, th0, toll)
% time_law Computes the true anomaly at a specific time in a PARABOLIC orbit, using 
% an inverse time law, given periapsis distance, gravitational parameter, and tolerance.
%
% PROTOTYPE
% [th, D, M] = time_law(t, rp, mu, t0, th0, toll)
%
% INPUT:
% t[nx1]         Array of time steps for which the true anomaly is calculated [T]
% rp[1]          Periapsis distance of the orbit [L]
% mu[1]          Gravitational parameter of the central body [L^3/T^2]
% t0[1]          Initial time [T]
% th0[1]         Initial true anomaly at time t0 [rad]
% toll[1]        Tolerance for the solution accuracy [-]
%
% OUTPUT:
% th[nx1]        Array of true anomaly values corresponding to each time step [rad]
% D[nx1]         Array of intermediate variable values for each time step [-]
% M[nx1]         Array of mean anomaly values for each time step [-]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-27: First version
%
% -------------------------------------------------------------------------

%% Initial Setup and Parameters
% Define initial values and orbital parameters used in calculations.

% Calculate initial D based on the initial true anomaly th0.
D0 = tan(th0 / 2);

% Compute semi-latus rectum, p, and specific angular momentum, h.
p = 2 * rp;
h = sqrt(mu * p);

% Calculate mean motion, n, based on the gravitational parameter and h.
n = mu^2 / h^3;

%% Calculate Mean Anomaly M
% Determine the mean anomaly M from the initial and target times.

% Mean anomaly M at each time step, adjusted by initial time t0.
M = n .* (t - t0);

%% Solve for Intermediate Variable D Using fsolve
% Use an optimization approach to solve for D by minimizing the function to meet the tolerance.

% Set options for fsolve to control tolerance and suppress display.
options = optimoptions('fsolve', 'FunctionTolerance', toll, 'Display', 'none');

% Solve for D using fsolve, finding the value of D that matches the mean anomaly M.
D = fsolve(@(D) 1/2 * (D + D.^3 / 3) - 1/2 * (D0 + D0^3 / 3) - M, 0, options);

%% Calculate True Anomaly th
% Convert the solved intermediate variable D back to the true anomaly th.

% True anomaly is twice the arctangent of D.
th = 2 * atan(D);

end
