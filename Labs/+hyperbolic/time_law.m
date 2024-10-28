function [th, F, M] = time_law(t, e, a, mu, t0, th0, toll)
% time_law Computes the true anomaly at a specific time in a HYPERBOLIC orbit
% based on the given time, eccentricity, semi-major axis, and gravitational parameter.
%
% PROTOTYPE
% [th, F, M] = time_law(t, e, a, mu, t0, th0, toll)
%
% INPUT:
% t[nx1]         Array of time values for which to compute the true anomaly [T]
% e[1]           Orbital eccentricity (e > 1 for hyperbolic orbits) [-]
% a[1]           Semi-major axis of the orbit [L]
% mu[1]          Gravitational parameter of the central body [L^3/T^2]
% t0[1]          Initial time [T]
% th0[1]         Initial true anomaly [rad]
% toll[1]        Tolerance for convergence in the fsolve function [-]
%
% OUTPUT:
% th[nx1]        Array of true anomaly values corresponding to each time step [rad]
% F[nx1]         Array of hyperbolic eccentric anomaly values [-]
% M[nx1]         Array of mean anomaly values for each time step [-]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-27: First version
%
% -------------------------------------------------------------------------

%% Initial Hyperbolic Eccentric Anomaly (F0)
% Compute the initial hyperbolic eccentric anomaly F0 from the given initial true anomaly th0.

F0 = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(th0 / 2));

%% Orbital Parameters Calculation
% Calculate orbital parameters, including the semi-latus rectum and mean motion.

% Semi-latus rectum of the hyperbolic trajectory.
p = a * (1 - e^2);

% Specific angular momentum.
h = sqrt(p * mu);

% Hyperbolic mean motion.
n = mu^2 / h^3 * (e^2 - 1)^(3/2);

%% Mean Anomaly (M) Calculation
% Calculate mean anomaly M as a function of time.

M = n * (t - t0);

%% Solve for Hyperbolic Eccentric Anomaly (F) Using fsolve
% Find F values that satisfy the mean anomaly equation using an iterative solver.

% Options for the fsolve function, set to desired tolerance.
options = optimoptions('fsolve', 'FunctionTolerance', toll, 'Display', 'none');

% Use fsolve to solve for F with initial guess of 0.
F = fsolve(@(F) e * sinh(F) - F - e * sinh(F0) + F0 - M, 0, options);

%% Compute True Anomaly (th) from F
% Convert the hyperbolic eccentric anomaly F back to the true anomaly th.

th = 2 * atan(sqrt((e + 1) / (e - 1)) * tanh(F / 2));

end
