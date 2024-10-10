function [alpha, delta, lon, lat, t] = groundTrack_mine(car, kep, Green_lon0, om_E, t, mu)
% groundTrack_mine - Computes the ground track of a spacecraft based on its orbital elements or Cartesian coordinates.
%
% PROTOTYPE:
% [alpha, delta, lon, lat, t] = groundTrack_mine(car, kep, Green_lon0, om_E, t, mu)
%
% INPUT:
% car [nx6]    - Cartesian state vector at each time step [rx, ry, rz, vx, vy, vz] [L, L/T]
% kep [1x6]    - Orbital elements [a, e, i, OM, w, th0] [L, -, rad]
% Green_lon0 [1] - Initial Greenwich longitude [rad]
% om_E [1]     - Angular velocity of Earth rotation [rad/s]
% t [nx1]      - Time vector [T]
% mu [1]       - Gravitational parameter of the primary body [L^3/T^2]
%
% OUTPUT:
% alpha [nx1]  - Right ascension [rad]
% delta [nx1]  - Declination [rad]
% lon [nx1]    - Longitude [rad]
% lat [nx1]    - Latitude [rad]
% t [nx1]      - Time vector with NaNs for discontinuities
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

    % Check if time vector is a row and transpose if necessary
    if size(t, 2) > 1
        t = t'; % Ensure t is a column vector
    end
    
    % Convert orbital elements or Cartesian coordinates to ECI coordinates
    if isempty(car)
        [r, v, t] = kep2ECI(t, kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu);
    elseif isempty(kep)
        if size(car, 2) > 6
            car = car'; % Transpose car to ensure correct dimensions
        end
        [r, v, t] = car2ECI(t, car(:, 1:3), car(:, 4:6), mu);
    end
    
    % Initial state data
    v_0 = v(1, :);  % Initial velocity vector
    r_0 = r(1, :);  % Initial position vector
    
    % Calculate semi-major axis (a) and orbital period (T)
    a = 1 / (2 / norm(r_0) - dot(v_0, v_0) / mu); % Semi-major axis
    T = 2 * pi * sqrt(a^3 / mu); % Orbital period

    % Calculate the Earth's rotation angle for each time instance
    theta_E = om_E * t; % Earth rotation during time t
    
    % Convert ECI coordinates to right ascension and declination
    [alpha, delta] = eci2ra(r);

    % Calculate longitude considering Earth's rotation
    lon = wrapTo2Pi(alpha - theta_E - Green_lon0); % Corrected longitude
    lat = delta; % Latitude remains unchanged

    % Detect large jumps in longitude to avoid discontinuities in the plot
    lon_diff = diff(lon); % Successive differences of longitude
    large_jump_lon = [false; abs(lon_diff) > pi]; % Detect jumps > 180 degrees

    % Insert NaN at discontinuity points to break the plot
    lon(large_jump_lon) = NaN;
    lat(large_jump_lon) = NaN;
    t(large_jump_lon) = NaN;

    % Add NaN at the end to close the plot correctly
    lon(end + 1) = NaN; 
    lat(end + 1) = NaN;
    t(end + 1) = NaN;
     
    % Load an image to use as the background (e.g., image of the Earth)
    img = imread('8081_earthmap10k.jpg');  % Ensure the file exists in the correct path
    
    % Display the image in the plot
    figure;
    imagesc([-180 180], [-90 90], flipud(img));  % Adjust image to geographic coordinates
    set(gca, 'YDir', 'normal');  % Correct Y-axis orientation
    
    % Ensure the image is not erased and overlay the plot
    hold on;
    
    % Plot the trajectory with a comet effect to visualize movement
    plot(rad2deg(lon) - 180, rad2deg(lat), 'Color', 'r', 'LineWidth', 2);
    plot(rad2deg(lon(1)) - 180, rad2deg(lat(1)), 'Color', 'b', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 13);
    plot(rad2deg(lon(end - 1)) - 180, rad2deg(lat(end - 1)), 'Color', 'g', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 13);
    xlim([-180 180]);
    ylim([-90 90]);
    xlabel('Longitude [deg]'); xticks(linspace(-180,180,13));
    ylabel('Latitude [deg]'); yticks(linspace(-90,90,7))
    hold off;
    grid on; ax = gca;
    ax.GridColor = 'k';
    
    % Create a comet plot for the trajectory
    %figure;
    %comet(rad2deg(lon) - 180, rad2deg(lat));
    % Set axis limits
    %xlim([-180 180]);
    %ylim([-90 90]);
    %hold off;
    %grid on;
end
