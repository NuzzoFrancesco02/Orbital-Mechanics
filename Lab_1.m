%% EX 1
%{
Orbit Propagation Code:
- Identify the states of the system and the relevant physical parameters.
- Write the second-order ODE to describe the dynamics (acceleration due to gravity).
- Reduce the second-order ODE into a first-order system (velocity and position).
- Implement the ODE function (odefun) for this system.
- Write the main script to numerically integrate the system using one of MATLAB's solvers,
  setting solver options as needed.

The equation of motion: r_dot_dot + mu/(norm(r)^3) * r = 0
where r_dot_dot is the acceleration, mu is the gravitational parameter, and r is the position.
%}

% Define the gravitational parameter for Earth (mu) [L^3/T^2].
mu = astroConstants(13);

% Initial conditions:
% r_0: Initial position vector (in km).
% v_0: Initial velocity vector (in km/s).
r_0 = [26578.137, 0, 0]'; 
v_0 = [0, 2.221, 3.173]';

% If needed, other initial conditions can be used (commented out below):
r_0 = [6495, -970, -3622]'; 
v_0 = [4.752, 2.130, 7.950]';

% Combine position and velocity into one state vector (s_0).
s_0 = [r_0; v_0];

% Calculate the norm (magnitude) of the initial velocity and position vectors.
v0_norm = norm(v_0); 
r0_norm = norm(r_0);

% Calculate the semi-major axis (a) based on orbital energy considerations.
a = 1 / (2 / norm(r_0) - dot(v_0, v_0) / mu);

% Compute the orbital period (T) using Kepler's third law.
T = 2 * pi * sqrt(a^3 / mu);

%% Non-perturbed Orbit Propagation
% Call the function to integrate the two-body problem without perturbations.
keplerian_orbit_ode(s_0, [0 5*T], mu);

%% Perturbed Orbit Propagation (with J2 effect)
% Include the J2 gravitational effect in the simulation:
% Re: Earth's mean radius [km].
% J2: Earth's J2 zonal harmonic (flattening effect).
Re = astroConstants(23); % Earth's radius [km]
J2 = astroConstants(9);  % J2 constant (unitless)

% Propagate the orbit with the J2 effect included.
keplerian_orbit_ode(s_0, [0 5*T], mu, 'J2', J2, Re);

%% Long-term Perturbed Orbit Propagation (600 orbits / ~1 year)
% Set numerical integration options with high precision.
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Number of orbital periods to display.
N_periods = 600;

% Time span for the integration (600 periods, with 1000 points per period).
tspan = linspace(0, N*T, N*1000);

% Numerically integrate the system including the J2 effect using ode113.
[t, y] = ode113(@(t, y) ode_2bp_J2(t, y, mu, J2, Re), tspan, s_0, options);

% Create a color map to change the orbit color over each period.
cmap = jet(N_periods);

% Plot the Earth as a 3D sphere.
Earth3d; 
grid on;

% Loop through each orbital period and plot the corresponding trajectory with a unique color.
for i = 1:N_periods
    % Indices that correspond to the current period.
    idx_start = (i-1) * 1000 + 1;
    idx_end = i * 1000;
    
    % Plot the segment of the orbit corresponding to the current period with its unique color.
    plot3(y(idx_start:idx_end, 1), y(idx_start:idx_end, 2), y(idx_start:idx_end, 3), ...
          'Color', cmap(i, :), 'LineWidth', 1); 
    hold on;
end
