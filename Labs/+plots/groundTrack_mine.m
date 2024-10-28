function [alpha, delta, lon, lat, t] = groundTrack_mine(t0, a, e, i, OM, om, th0, thf, Green_lon0, om_E, mu, str)
% groundTrack_mine - Computes the ground track of a spacecraft based on orbital elements.
%
% PROTOTYPE:
% [alpha, delta, lon, lat, t] = groundTrack_mine(t0, a, e, i, OM, om, th0, thf, Green_lon0, om_E, mu, str)
%
% INPUT:
% t0           - Initial time [T]
% a            - Semi-major axis [L]
% e            - Eccentricity [-]
% i            - Inclination [rad]
% OM           - RAAN (Right Ascension of Ascending Node) [rad]
% om           - Argument of periapsis [rad]
% th0          - Initial true anomaly [rad]
% thf          - Final true anomaly [rad]
% Green_lon0   - Initial Greenwich longitude [rad]
% om_E         - Angular velocity of Earth rotation [rad/s]
% mu           - Gravitational parameter [L^3/T^2]
% str          - Optional parameter to consider J2 perturbation
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
% -------------------------------------------------------------------------

    % Generate the true anomaly vector and compute the time vector
    th = linspace(th0, thf, 10000);
    t = elliptic.kepler_inv_equation(e, a, mu, t0, th);
    
    % Apply J2 perturbation if selected
    J2 = astroConstants(9);
    Re = astroConstants(23);
    if strcmp(str, 'J2')
        fac = -3/2 * sqrt(mu) * J2 * Re^2 / (1 - e^2)^2 / a^(7/2);
    else
        fac = 0;
    end
    OMdot = fac * cos(i);
    omdot = fac * (5/2 * sin(i)^2 - 2);

    % Compute position in ECI, transform to ECEF
    r_rel = [];
    for j = 1 : length(t)
        r = coord.kep2car_theta(a, e, i, OM + OMdot * t(j), om + omdot * t(j), th(j), mu);
        theta = om_E * (t(j) - t(1));
        Q = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
        r_rel = [r_rel; (Q * r)'];
    end

    % Convert to right ascension and declination
    [alpha, delta] = coord.eci2ra(r_rel);

    % Calculate longitude and latitude
    lon = wrapTo2Pi(alpha - Green_lon0);
    lat = delta;

    % Handle large jumps in longitude
    lon_diff = diff(lon);
    large_jump_lon = [false; abs(lon_diff) > pi / 2];
    lon(large_jump_lon) = NaN;
    lat(large_jump_lon) = NaN;
    t(large_jump_lon) = NaN;

    % Display the ground track over an Earth map
    img = imread('8081_earthmap10k.jpg');
    figure;
    imagesc([-180 180], [-90 90], flipud(img));
    set(gca, 'YDir', 'normal');
    hold on;
    
    % Plot the trajectory
    plot(rad2deg(lon) - 180, rad2deg(lat), 'r', 'LineWidth', 2);
    plot(rad2deg(lon(1)) - 180, rad2deg(lat(1)), 'bo', 'MarkerSize', 13, 'LineWidth', 2);
    plot(rad2deg(lon(end - 1)) - 180, rad2deg(lat(end - 1)), 'go', 'MarkerSize', 13, 'LineWidth', 2);
    xlim([-180 180]);
    ylim([-90 90]);
    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    grid on;
    hold off;
    
    % Comet plot for the trajectory
    %figure;
    % comet(rad2deg(lon) - 180, rad2deg(lat));
    % xlim([-180 180]);
    % ylim([-90 90]);
    % grid on;
end
