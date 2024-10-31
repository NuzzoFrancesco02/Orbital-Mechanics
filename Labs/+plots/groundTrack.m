function [alpha, delta, lon, lat, t] = groundTrack(t, r, v, n_orb, Green_lon0, om_E, mu, str)
% groundTrack_mine - Computes the ground track of a spacecraft based on orbital elements.
% This function calculates and visualizes the ground track of a spacecraft,
% transforming its position from ECI (Earth-Centered Inertial) coordinates 
% to latitudes and longitudes on the Earth's surface. Optional J2 perturbations 
% can be applied to adjust the orbital elements over time.
%
% PROTOTYPE:
% [alpha, delta, lon, lat, t] = groundTrack_mine(t0, r, v, n_orb, Green_lon0, om_E, mu, str)
%
% INPUT:
% t            - Initial time or time array when r is an array [T]
% r            - Initial position vector in ECI or collection of position in ECI [L]
% v            - Initial velocity vector in ECI [L/T]
% n_orb        - Number of orbits to simulate [-]
% Green_lon0   - Initial Greenwich longitude [rad]
% om_E         - Earth's rotation rate [rad/s]
% mu           - Gravitational parameter of the central body [L^3/T^2]
% str          - Optional parameter ('J2') to consider J2 perturbation (default: no J2 effect)
%
% OUTPUT:
% alpha        - Right ascension [rad]
% delta        - Declination [rad]
% lon          - Longitude [rad]
% lat          - Latitude [rad]
% t            - Time vector with NaNs for discontinuities
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-28: Initial version
%
% -------------------------------------------------------------------------
    t0 = t(1);
    % If initial input is in Cartesian form, convert to Keplerian elements
    if min(size(r)) == 1 && min(size(v)) == 1
        [a, e, i, OM, om, th0] = coord.car2kep_theta(r, v, mu); % Cartesian to Keplerian
        T = 2 * pi * sqrt(a^3 / mu);                           % Calculate orbit period
        thf = elliptic.kepler_equation(n_orb * T, e, a, mu, t0, th0, 1e-14); % Final true anomaly
        th = linspace(th0, thf, 20000);                         % True anomaly vector
        t = elliptic.kepler_inv_equation(e, a, mu, t0, th);    % Time vector from anomaly
        J2 = astroConstants(9);                                % J2 perturbation constant
        Re = astroConstants(23);                               % Earth's radius
        if strcmp(str, 'J2')
            fac = -3/2 * sqrt(mu) * J2 * Re^2 / (1 - e^2)^2 / a^(7/2);
        else
            fac = 0;
        end
        OMdot = fac * cos(i);                                  % RAAN rate with J2 effect
        omdot = fac * (5/2 * sin(i)^2 - 2);                    % Perigee rate with J2 effect
    else
        r_rel = r;
    end

    % Compute ECI to ECEF transformation and get relative position vector
    if min(size(r)) == 1 && min(size(v)) == 1
        r_rel = [];
        for j = 1 : length(t)
            r = coord.kep2car_theta(a, e, i, OM + OMdot * t(j), om + omdot * t(j), th(j), mu);
            theta = om_E * (t(j) - t(1));                     % Earth rotation angle
            Q = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1]; % Rotation matrix
            r_rel = [r_rel; (Q * r)'];
        end
    end

    % Convert relative position to right ascension and declination
    [alpha, delta] = coord.eci2ra(r_rel);
    
    if min(size(r)) ~= 1 && min(size(v)) ~= 1
        if size(t,2) > 3
            t = t';
        end
        alpha = alpha + om_E .* (t-t(1));
    end

    % Calculate geographic longitude and latitude
    lon = wrapToPi(alpha - Green_lon0);
    lat = delta;

    % Handle large jumps in longitude for continuity in plotting
    lon_diff = diff(lon);
    large_jump_lon = [false; abs(lon_diff) > pi / 2];
    lon(large_jump_lon) = NaN;
    lat(large_jump_lon) = NaN;
    t(large_jump_lon) = NaN;

       % Plotting the ground track on a world map background
    h = figure('Name', 'Ground Track');
    set(h, 'Units', 'Normalized', 'OuterPosition', [.15 .25 .7 .7]);
    axis equal;
    set(gca, 'XTick', -180:15:180, 'YTick', -90:10:90, 'XTickMode', 'manual', 'YTickMode', 'manual');
    
    % Display Earth map image as background
    image_file = '8081_earthmap10k.jpg';
    cdata = flip(imread(image_file));
    h_img = imagesc([-180, 180], [-90, 90], cdata);
    set(gca, 'YDir', 'normal'); % Imposta la direzione dell'asse y normale
    hold on;
    
    xlim([-180, 180]);
    ylim([-90, 90]);
    
    set(h_img, 'AlphaData', 0.9); % Adjust transparency for visibility
    plot(rad2deg(lon), rad2deg(lat), 'r', 'LineWidth', 1); % Ground track in red
    plot(rad2deg(lon(1)), rad2deg(lat(1)), 'go', 'MarkerSize', 13, 'LineWidth', 2); % Start point in green
    plot(rad2deg(lon(end - 1)), rad2deg(lat(end - 1)), 'g', 'MarkerSize', 13, 'LineWidth', 2, 'Marker', 'square'); % End point
    xlabel('Longitude \lambda [deg]');
    ylabel('Latitude \phi [deg]');
    legend({'Ground Track', 'Start Position', 'End Position'}, 'Location', 'best');
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridLineWidth', 1, 'GridColor', 'k');
    hold off;
end