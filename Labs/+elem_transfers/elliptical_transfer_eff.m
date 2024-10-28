function elliptical_transfer_eff(e1, mu)
    % elliptical_transfer_eff Plots the delta-v efficiency ratio for an elliptical orbit transfer.
    %
    % PROTOTYPE
    % elliptical_transfer_eff(a1, e1, mu)
    %
    % INPUT:
    % a1[1]        Semi-major axis of the initial orbit [L]
    % e1[1]        Eccentricity of the initial orbit [-]
    % mu[1]        Gravitational parameter of the primary body [L^3/T^2]
    %
    % OUTPUT:
    % None (function generates plots)
    %
    % CONTRIBUTORS:
    % Francesco Nuzzo
    %
    % VERSIONS:
    % 2024-10-27: First version
    %
    % -------------------------------------------------------------------------

    %% Orbital Parameters
    % Calculating the apogee and perigee of the initial orbit.
    a1 = 1e6;
    rp = a1 * (1 - e1); % Perigee distance in the initial orbit
    ra = a1 * (1 + e1); % Apogee distance in the initial orbit

    % Set rA and rA_prime to represent the initial and final apogee distances.
    rA = rp;            % Redefine rA as the perigee distance of the initial orbit
    rA_prime = ra;      % Redefine rA_prime as the apogee distance of the initial orbit

    %% Delta-V Calculation Functions
    % Define functions to calculate delta-v for the transfer based on input parameters alpha and beta.
    D_v = @(alpha, beta) abs(sqrt(2 * mu) ./ rA .* (-sqrt(rA * rA_prime ./ (rA + rA_prime)) + sqrt(rA ./ (1 ./ alpha + 1)))) + ...
                         abs(sqrt(2 * mu) ./ (alpha * rA) .* (sqrt((rA .* alpha) ./ (1 + alpha ./ beta)) - sqrt(rA ./ (1 ./ alpha + 1))));

    gamma = (e1 + 1) / (1 - e1); % Compute the eccentricity ratio adjustment.
    D_v_prime = @(alpha, beta) abs(sqrt(2 * mu) ./ rA_prime .* (-sqrt(rA .* rA_prime ./ (rA + rA_prime)) + sqrt(rA_prime ./ (gamma ./ beta + 1)))) + ...
                               abs(sqrt(2 * mu) ./ (beta * rA) .* (sqrt((rA .* alpha) ./ (1 + alpha ./ beta)) - sqrt(rA_prime ./ (gamma ./ beta + 1))));

    D_v_rapp = @(alpha, beta) D_v_prime(alpha, beta) ./ D_v(alpha, beta); % Delta-v ratio

    %% Grid Generation for Plotting

    [X, Y] = meshgrid(linspace(1, 10, 100), linspace(1, 10, 100)); % Generate grid values

    Z = D_v_rapp(X, Y); % Calculate delta-v ratio for efficiency

    %% Plot Results
    % Surface plot of delta-v efficiency ratio.
    figure;
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.7); % Use FaceAlpha for better visibility
    hold on;

    % Add a plane at Z = 1
    Z_plane = ones(size(X)); % Create a matrix for the Z = 1 plane
    surf(X, Y, Z_plane, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'k'); % Plot the Z=1 plane

    colorbar; % Add color bar for reference
    colormap(pink); % Change color map for better visualization
    grid on;
    title('Delta-V Efficiency Ratio (Surface Plot)', 'FontWeight', 'bold', 'FontSize', 14);
    xlabel('\alpha', 'FontWeight', 'bold');
    ylabel('\beta', 'FontWeight', 'bold');
    zlabel('D\_v\_rapp(\alpha, \beta)', 'FontWeight', 'bold');
    axis tight; % Tighten axis limits for better focus
    set(gca, 'FontSize', 12); % Set font size for axes

    % Contour plot of delta-v efficiency ratio.
    figure;
    % Include a contour line for the value of 1 explicitly
    [C, h] = contour(X, Y, Z, 'LineWidth', 1.5); % General contour lines
    hold on; % Hold the plot for additional contours
    contour(X, Y, Z, [1 1], 'LineColor', 'k', 'LineWidth', 2); % Explicit contour line for Z = 1
    clabel(C, h, 'FontSize', 10); % Label the contours
    colorbar; % Add color bar for reference
    colormap(pink); % Change color map for contour
    grid on;
    title('Delta-V Efficiency Ratio (Contour Plot)', 'FontWeight', 'bold', 'FontSize', 14);
    xlabel('\alpha', 'FontWeight', 'bold');
    ylabel('\beta', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12); % Set font size for axes
    axis tight; % Tighten axis limits for better focus

    hold off; % Release the hold
end
