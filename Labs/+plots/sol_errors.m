function [err_h, err_e, err_perp, err_eps] = sol_errors(t, r, v, mu)
% This function computes the errors in angular momentum, eccentricity, 
% perpendicularity condition, and specific energy over time, and plots these errors.
%
% PROTOTYPE:
% [err_h, err_e, err_perp, err_eps] = sol_errors(t, r, v, mu)
%
% INPUT:
% t [nx1]    - Time intervals [T]
% r [nx3]    - Position vector at each time step (rx, ry, rz) [L]
% v [nx3]    - Velocity vector at each time step (vx, vy, vz) [L/T]
% mu [1]     - Gravitational parameter of the primary body [L^3/T^2]
%
% OUTPUT:
% err_h [nx1]    - Error in angular momentum: ||h(t) - h(0)|| [L^2/T]
% err_e [nx1]    - Error in eccentricity vector: ||e(t) - e(0)|| [-]
% err_perp [nx1] - Perpendicularity condition error: ||dot(h(t), e(t))|| [-]
% err_eps [nx1]  - Specific orbital energy error: |eps(t) - eps(0)| [L^2/T^2]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

    % Ensure that the position and velocity vectors are column-major
    [row, col] = size(r);
    if row < col
        r = r';  % Transpose if needed
    end
    [row, col] = size(v);
    if row < col
        v = v';  % Transpose if needed
    end

    % Compute angular momentum vector at each time step
    h = cross(r, v, 2);  % Cross product between position and velocity
    h_0 = h(1,:);        % Initial angular momentum
    err_h = vecnorm(h - h_0, 2, 2);  % Error in angular momentum magnitude

    % Compute eccentricity vector at each time step
    e = cross(v, cross(r, v, 2), 2) ./ mu - r ./ vecnorm(r, 2, 2);  % Orbital eccentricity
    e_0 = e(1,:);        % Initial eccentricity
    err_e = vecnorm(e - e_0, 2, 2);  % Error in eccentricity magnitude

    % Perpendicularity condition: dot product of h and e should be zero
    err_perp = dot(h, e, 2);  % Error in perpendicularity condition

    % Compute specific orbital energy at each time step
    eps_t = vecnorm(v, 2, 2).^2 ./ 2 - mu ./ vecnorm(r, 2, 2);  % Specific orbital energy
    eps_0 = eps_t(1,:);   % Initial specific orbital energy
    err_eps = vecnorm(eps_t - eps_0, 2, 2);  % Error in specific orbital energy

    % Plot the errors in a 2x2 grid of subplots
    figure()

    % Plot angular momentum error
    subplot(2, 2, 1)
    plot(t, err_h, 'LineWidth', 2); 
    title('h_{vec} error'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('||h(t)-h_0|| [km^2/s]');
    fprintf('Max error for ||h(t)-h_0|| = %.4e\n', max(err_h));
    xlim([t(1) t(end)])

    % Plot eccentricity error
    subplot(2, 2, 2)
    plot(t, err_e, 'LineWidth', 2); 
    title('e_{vec} error'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('||e(t)-e_0|| [-]');
    fprintf('Max error for ||e(t)-e_0|| = %.4e\n', max(err_e));
    xlim([t(1) t(end)])

    % Plot perpendicularity condition error
    subplot(2, 2, 3)
    plot(t, err_perp, 'LineWidth', 2); 
    title('perp. condition error'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('||dot(h(t),e(t))|| [km^2/s]');
    fprintf('Max error for ||dot(h(t),e(t))|| = %.4e\n', max(err_perp));
    xlim([t(1) t(end)])

    % Plot specific energy error
    subplot(2, 2, 4)
    plot(t, err_eps, 'LineWidth', 2); 
    title('energy error'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('|eps(t)-eps_0| [km^2/s^2]');
    fprintf('Max error for |eps(t)-eps_0| = %.4e\n', max(err_eps));
    xlim([t(1) t(end)])

    % Plot the errors in a 2x2 grid of subplots
    figure()

    % Plot angular momentum error
    subplot(2, 2, 1)
    plot(t, h, 'LineWidth', 2); 
    title('h_{vec}'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('h(t) [km^2/s]');
    xlim([t(1) t(end)])

    % Plot eccentricity error
    subplot(2, 2, 2)
    plot(t, e, 'LineWidth', 2); 
    title('e_{vec}'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('e(t) [-]');
    xlim([t(1) t(end)])

    % Plot perpendicularity condition error
    subplot(2, 2, 3)
    plot(t, err_perp, 'LineWidth', 2); 
    title('perp. condition error'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('||dot(h(t),e(t))|| [km^2/s]');
    fprintf('Max error for ||dot(h(t),e(t))|| = %.4e\n', max(err_perp));
    xlim([t(1) t(end)])

    % Plot specific energy error
    subplot(2, 2, 4)
    plot(t, eps_t, 'LineWidth', 2); 
    title('energy'); 
    grid on; 
    xlabel('t [s]'); 
    ylabel('eps(t) [km^2/s^2]');
    xlim([t(1) t(end)])

end
