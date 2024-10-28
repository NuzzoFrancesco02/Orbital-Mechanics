function [v_r, v_t, u_r, u_h, u_t] = velocity_components(t,r,v)
% This function computes the radial and transversal components of velocity
% and plots them over time.
%
% PROTOTYPE:
% [v_r, v_t, u_r, u_h, u_t] = velocity_components(t, r, v)
%
% INPUT:
% t [nx1]    - Time intervals [T]
% r [nx3]    - Position vector at each time step (rx, ry, rz) [L]
% v [nx3]    - Velocity vector at each time step (vx, vy, vz) [L/T]
%
% OUTPUT:
% v_r        - Radial component of velocity at each time step [L/T]
% v_t        - Transversal component of velocity at each time step [L/T]
% u_r        - Radial unit vector at each time step [-]
% u_h        - Unit vector parallel to angular momentum at each time step [-]
% u_t        - Transversal unit vector at each time step [-]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

% Compute the angular momentum vector at each time step using the cross product
% of position and velocity vectors
h = cross(r, v, 2);

% Compute the radial unit vector (u_r) by normalizing the position vector r
u_r = r ./ vecnorm(r, 2, 2);

% Compute the unit vector u_h parallel to the angular momentum vector h
u_h = h ./ vecnorm(h, 2, 2);

% Compute the transversal unit vector u_t by taking the cross product
% of u_h and u_r
u_t = cross(u_h, u_r, 2);

% Compute the radial velocity component (v_r) by taking the dot product
% of the velocity vector and the radial unit vector u_r
v_r = dot(v, u_r, 2);

% Compute the transversal velocity component (v_t) by taking the dot product
% of the velocity vector and the transversal unit vector u_t
v_t = dot(v, u_t, 2);

% Plot the radial and transversal velocity components over time

% Create a new figure
figure()

% Plot the radial velocity component
subplot(2, 1, 1)
plot(t, v_r, 'LineWidth', 2); 
grid on; 
title('Radial Velocity v_{rad}');
xlim([t(1) t(end)]);

% Plot the transversal velocity component
subplot(2, 1, 2)
plot(t, v_t, 'LineWidth', 2); 
grid on; 
title('Transversal Velocity v_{trans}');
xlim([t(1) t(end)]);
