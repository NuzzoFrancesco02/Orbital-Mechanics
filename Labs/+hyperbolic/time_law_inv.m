function [t, F, M] = time_law_inv(e, a, mu, t0, th)
% time_law_inv Computes the time at a specific true anomaly in a HYPERBOLIC orbit
% using the inverse time law for hyperbolic trajectories.
%
% PROTOTYPE
% [t, F, M] = time_law_inv(e, a, mu, t0, th)
%
% INPUT:
% e[1]           Orbital eccentricity (e > 1 for hyperbolic orbits) [-]
% a[1]           Semi-major axis of the orbit [L]
% mu[1]          Gravitational parameter of the central body [L^3/T^2]
% t0[1]          Initial time [T]
% th[nx1]        Array of true anomaly values [rad]
%
% OUTPUT:
% t[nx1]         Array of time values corresponding to each true anomaly [T]
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

%% Orbital Parameters Calculation
% Compute orbital parameters based on the given eccentricity, semi-major axis, and gravitational parameter.

% Semi-latus rectum of the hyperbolic trajectory.
p = a * (1 - e^2);

% Specific angular momentum.
h = sqrt(p * mu);

% Hyperbolic mean motion.
n = mu^2 / h^3 * (e^2 - 1)^(3/2);

%% Hyperbolic Eccentric Anomaly F Calculation
% Compute F, the hyperbolic eccentric anomaly, from the true anomaly th.

% Conversion from true anomaly th to hyperbolic eccentric anomaly F.
F = 2 * atanh(sqrt((e - 1) / (e + 1)) .* tan(th / 2));

%% Mean Anomaly M Calculation
% Calculate the mean anomaly M from F and the eccentricity.

% Initial mean anomaly M0 at the first true anomaly in th.
M0 = e * sinh(F(1)) - F(1);

% Mean anomaly M at each true anomaly in the input array.
M = e * sinh(F) - F;

%% Time Calculation
% Calculate the time for each true anomaly based on M and the mean motion n.

% Calculate time t based on mean anomaly change relative to initial time t0.
t = t0 + (M - M0) / n;

end
