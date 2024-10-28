function [t, D, M] = time_law_inv(rp, mu, t0, th)
% time_law_inv Calculates the time evolution for a PARABOLIC orbit based on an 
% inverse time law, given the periapsis distance and gravitational parameter.
%
% PROTOTYPE
% [t, D, M] = time_law_inv(rp, mu, t0, th)
%
% INPUT:
% rp[1]        Periapsis distance of the orbit [L]
% mu[1]        Gravitational parameter of the central body [L^3/T^2]
% t0[1]        Initial time [T]
% th[nx1]      Array of true anomaly angles (in radians) for each time step [rad]
%
% OUTPUT:
% t[nx1]       Time array corresponding to each true anomaly angle [T]
% D[nx1]       Array of intermediate variable values for each true anomaly [-]
% M[nx1]       Array of mean anomaly values for each true anomaly angle [-]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-27: First version
%
% -------------------------------------------------------------------------

%% Parameter Calculation
% Define semi-latus rectum (p) and specific angular momentum (h) of the orbit.

% The semi-latus rectum, p, is calculated as twice the periapsis distance.
p = 2 * rp;

% Calculate the specific angular momentum, h, based on p and the gravitational parameter, mu.
h = sqrt(p * mu);

% Calculate mean motion, n, which is the rate at which the orbit progresses in terms of mean anomaly.
n = mu^2 / h^3;

%% Calculate D and Mean Anomaly M
% Use the true anomaly to calculate intermediate variable D and mean anomaly M.

% Intermediate variable D, which is tangent of half the true anomaly angle.
D = tan(th ./ 2);

% Mean anomaly M is calculated using D, representing the position along the orbit in terms of time.
M = 1/2 * (D + D.^3 / 3);

% Initial mean anomaly, M0, is the first value in M, used to calculate time evolution.
M0 = M(1);

%% Time Calculation
% Calculate the time corresponding to each true anomaly based on M and the initial time.

% Calculate time t by adjusting the mean anomaly with the initial mean anomaly M0.
t = t0 + (M - M0) / n;

end

