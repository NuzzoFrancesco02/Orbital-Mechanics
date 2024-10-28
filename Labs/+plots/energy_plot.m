function energy_plot(r_0, v_0, mu)
% Plots the specific orbital energy as a function of the distance r.
%
% This function computes and plots the specific orbital energy (eps) 
% based on the initial position, velocity, and gravitational parameter. 
% It also provides information about the type of orbit (elliptical, 
% parabolic, or hyperbolic) based on the value of eps_0.
%
% PROTOTYPE:
% energy_plot(r_0, v_0, mu)
%
% INPUT:
% r_0 [nx3]     - Initial position vector [rx, ry, rz] [L]
% v_0 [nx3]     - Initial velocity vector [vx, vy, vz] [L/T]
% mu [1]        - Gravitational parameter of the primary body [L^3/T^2]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

    % Calculate the specific energy at the initial conditions
    eps_0 = norm(v_0)^2 / 2 - mu / norm(r_0);  % Specific energy [L^2/T^2]
    
    % Calculate the specific angular momentum
    h_0 = cross(r_0, v_0);  % Specific angular momentum vector [L^2/T]
    
    % Define a range of distances for plotting the energy
    r_energy = linspace(1e3, norm(r_0) * 10, 10000);  % Distance range [L]
    
    % Compute specific energy as a function of distance r
    eps_r = norm(h_0)^2 / 2 ./ r_energy.^2 - mu ./ r_energy;
    
    % Find the minimum energy value
    eps_min = min(eps_r);

    % Classify the type of orbit based on eps_0
    if eps_0 < 0
        fprintf('\nElliptical orbit: eps_0 = %.4f < 0, while eps_min = %.4f\n', eps_0, eps_min);
    elseif eps_0 == 0
        fprintf('\nParabolic orbit: eps_0 = 0\n');
    elseif eps_0 > 0
        fprintf('\nHyperbolic orbit: eps_0 = %.4f > 0\n', eps_0);
    elseif eps_0 == eps_min
        fprintf('\nCircular orbit: eps_0 = eps_min = %.4f\n', eps_0);
    end
    
    %% Plotting the specific energy as a function of distance r
    figure()
    plot(r_energy, eps_r, 'LineWidth', 2);  % Plot energy curve
    grid on;
    hold on;
    
    % Add a horizontal line at eps = 0 for reference (parabolic case)
    yline(0, 'LineWidth', 1, 'Color', 'k');
    
    % Set plot limits
    ylim([eps_min * 1.2, 10]);
    xlim([r_energy(1), r_energy(end)]);
    
    % Mark the initial position energy on the plot
    plot(norm(r_0), eps_0, 'Marker', 'x', 'MarkerSize', 10, 'LineWidth', 3);
    
    % Add a dashed line for the minimum energy value
    yline(eps_min, 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Add a legend to explain the curves
    legend('eps(r)', '', 'eps_0', 'eps_{min}');
    
end
