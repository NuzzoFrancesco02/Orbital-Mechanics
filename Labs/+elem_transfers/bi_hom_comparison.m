function bi_hom_comparison
    % bi_hom_comparison Compares the velocities of two orbital models.
    %
    % This function computes the difference between the velocity
    % models D_v_H and D_v_B over specified ranges of the parameters 
    % alpha and beta. It visualizes the results using both a surface 
    % plot and a contour plot, indicating where the difference is zero.
    %
    % PROTOTYPE:
    % bi_hom_comparison()
    %
    % INPUT:
    % None
    %
    % OUTPUT:
    % A figure showing the surface difference and a contour plot
    %
    % CONTRIBUTORS:
    % Francesco Nuzzo
    %
    % VERSIONS:
    % 2024-10-27: Updated version with detailed comments
    %
    % -------------------------------------------------------------------------

    % Define the functions for D_v_H and D_v_B
    D_v_H = @(alpha) 1 ./ sqrt(alpha) - sqrt(2) * (1 - alpha) ./ (sqrt(alpha .* (1 + alpha))) - 1;
    D_v_B = @(alpha, beta) sqrt(2 .* (alpha + beta) ./ (alpha .* beta)) - (1 + sqrt(alpha)) ./ sqrt(alpha) - sqrt(2 ./ (beta .* (1 + beta))) .* (1 - beta);
    
    % Create a grid of values for alpha and beta
    alpha = linspace(5, 40, 100); % Range for alpha from 5 to 40
    beta = linspace(5, 100, 100); % Range for beta from 5 to 100
    [X, Y] = meshgrid(alpha, beta); % Create a meshgrid for plotting

    % Calculate values for D_v_H and D_v_B
    Z_H = D_v_H(X); % Compute D_v_H for the grid
    Z_B = D_v_B(X, Y); % Compute D_v_B for the grid

    % Calculate the difference between the two models
    Z_diff = Z_H - Z_B;

    % Create the figure for the surface plot
    figure;
    
    % Plot the surface of the difference
    surf(X, Y, Z_diff, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    hold on;

    % Plot the plane where Z=0 for reference
    Z_plane = zeros(size(X)); % Create a plane at Z=0
    surf(X, Y, Z_plane, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'k'); % Plot the plane

    % Set up the graph labels and title
    xlabel('\alpha'); % Label for the x-axis
    ylabel('\beta'); % Label for the y-axis
    zlabel('Difference D_v_H - D_v_B'); % Label for the z-axis
    title('Difference between D_v_H and D_v_B with intersection at Z=0'); % Title of the graph
    grid on; % Enable grid
    view(3); % Set view to 3D
    colormap(turbo)
    colorbar; % Add a color bar for indicating values
    hold off;

    % Create a figure for the contour plot
    figure();
    
    % Plot the contour of the difference
    [C, h] = contour(X, Y, Z_diff, 'LineWidth', 1.5); % Draw contours
    clabel(C, h, 'FontSize', 10); % Label the contour levels
    colormap(turbo)
    colorbar; % Add a color bar
    xlabel('\alpha'); % Label for the x-axis
    ylabel('\beta'); % Label for the y-axis
    title('Contour of the Difference D_v_H - D_v_B'); % Title for the contour plot
    grid on; % Enable grid
end
