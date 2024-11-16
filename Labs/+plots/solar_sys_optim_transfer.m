function solar_sys_optim_transfer(planet1, planet2, t_dep, t_arr, fac)
% solar_sys_optim_transfer - Visualizes the transfer orbit between two planets in the solar system.
%
% PROTOTYPE:
% solar_sys_optim_transfer(planet1, planet2, t_dep, t_arr, fac)
%
% INPUT:
% planet1 [int]  - Identifier of the departure planet (1: Mercury, ..., 9: Pluto) [-]
% planet2 [int]  - Identifier of the arrival planet (1: Mercury, ..., 9: Pluto) [-]
% t_dep   [str]  - Departure date in datetime format [T]
% t_arr   [str]  - Arrival date in datetime format [T]
% fac     [1]    - Scaling factor for visualization purposes [-]
%
% OUTPUT:
% None. The function plots the transfer trajectory and planetary orbits in a 3D visualization.
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

% Define planet names for labels
strs = {'Mercury','Venus','Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'};

% Convert input dates to Modified Julian Date 2000 (MJD2000)
t1 = timeConversion.date2mjd2000(t_dep);
t2 = timeConversion.date2mjd2000(t_arr);

% Retrieve solar system center (Sun) data
[kep_sun, mu] = uplanet(t1, 10); % Planet 10 corresponds to the Sun
r_sun = coord.kep2car_theta(kep_sun(1), kep_sun(2), kep_sun(3), kep_sun(4), kep_sun(5), kep_sun(6), mu);

% Plot the Milky Way background and the Sun

background('Milky Way');
hold on;  % Keep everything in the plot

opts = struct('Position', r_sun, 'Units', 'km', 'RefPlane', 'ecliptic');
sun = planet3D('Sun', opts); drawnow;
sun.XData = sun.XData .* 1e-2; % Scale for visualization
sun.YData = sun.YData .* 1e-2; 
sun.ZData = sun.ZData .* 1e-2; 

% Plot departure planet orbit and position
[kep1, ~] = uplanet(t1, planet1);
kep2 = uplanet(t2, planet1);
if kep2(6) < kep1(6)
    kep2(6) = kep2(6)+2*pi;
end
r1 = coord.kep2car_theta(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), linspace(kep1(6), kep2(6), 1000), mu);
r1_complete = coord.kep2car_theta(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), linspace(0, 2*pi, 1000), mu);

% Plot departure planet trajectory
opts = struct('Clipping', 'on', 'Position', r1(1, :) .* fac, 'Units', 'km', 'RefPlane', 'equatorial'); 
planet3D(strs{planet1}, opts); drawnow;
opts = struct('Clipping', 'on', 'Position', r1(end, :) .* fac, 'Units', 'km', 'RefPlane', 'equatorial');
planet3D(strs{planet1}, opts); drawnow;

% Plot departure planet's orbit
h1 = plot3(r1(:, 1) .* fac, r1(:, 2) .* fac, r1(:, 3) .* fac, 'LineWidth', 2, 'Color', "#77AC30"); drawnow;
plot3(r1_complete(:, 1) .* fac, r1_complete(:, 2) .* fac, r1_complete(:, 3) .* fac, 'LineWidth', 1, 'LineStyle', '--', 'Color', "#77AC30"); drawnow;

% Add departure date label
text((r1(1, 1) + 2e3*astroConstants(planet1+20)) * fac, ...
    (r1(1, 2) + 2e3*astroConstants(planet1+20)) * fac, ...
    (r1(1, 3) + 2e3*astroConstants(planet1+20)) * fac, ...
    ['Dep: ', datestr(t_dep, 'yyyy-mm-dd')], ...
    'Color', "#77AC30", 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'none');

% Plot arrival planet orbit and position
[kep1, ~] = uplanet(t1, planet2);
kep2 = uplanet(t2, planet2);
if kep2(6) < kep1(6)
    kep2(6) = kep2(6)+2*pi;
end
r2 = coord.kep2car_theta(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), linspace(kep1(6), kep2(6), 1000), mu);
r2_complete = coord.kep2car_theta(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), linspace(0, 2*pi, 1000), mu);

opts = struct('Position', r2(1, :) .* fac, 'Units', 'km', 'RefPlane', 'equatorial');
planet3D(strs{planet2}, opts); drawnow;
opts = struct('Position', r2(end, :) .* fac, 'Units', 'km', 'RefPlane', 'equatorial');
planet3D(strs{planet2}, opts); drawnow;

% Plot arrival planet's orbit
h2 = plot3(r2(:, 1) .* fac, r2(:, 2) .* fac, r2(:, 3) .* fac, 'LineWidth', 2, 'Color', "#A2142F"); 
plot3(r2_complete(:, 1) .* fac, r2_complete(:, 2) .* fac, r2_complete(:, 3) .* fac, 'LineWidth', 1, 'LineStyle', '--', 'Color', "#A2142F"); drawnow;

% Add arrival date label
text((r2(end, 1) + 2e3*astroConstants(planet2+20))* fac, ...
    (r2(end, 2) + 2e3*astroConstants(planet2+20))* fac, ...
    (r2(end, 3) + 2e3*astroConstants(planet2+20)) * fac, ...
    ['Arr: ', datestr(t_arr, 'yyyy-mm-dd')], ...
    'Color', "#A2142F", 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'none');

% Compute Lambert transfer orbit
Nrev = 0; optionsLMR = 1; orbitType = 0; Ncase = 0;
[at, ~, et, ~, vt1, vt2, ~, ~] = lambert.lambertMR(r1(1, :), r2(end, :), (t2 - t1) * 24 * 60 * 60, mu, orbitType, Nrev, Ncase, optionsLMR);

% Convert initial and final transfer orbit states to Keplerian elements
[~, ~, it, OMt, omt, th1t] = coord.car2kep_theta(r1(1, :), vt1, mu);
[~, ~, ~, ~, ~, th2t] = coord.car2kep_theta(r2(end, :), vt2, mu);

% Plot transfer orbit
if th2t < th1t
    th2t = th2t + 2*pi;
end
rt = coord.kep2car_theta(at, et, it, OMt, omt, linspace(th1t, th2t, 100), mu);
h3 = plot3(rt(:, 1) .* fac, rt(:, 2) .* fac, rt(:, 3) .* fac, 'LineWidth', 2, 'Color', "#EDB120"); drawnow;

% Add legends and labels
legend([h1, h2, h3], {strcat(strs{planet1}, ' Orbit'), strcat(strs{planet2}, ' Orbit'), 'Transfer Orbit'}, 'Location', 'best');

title(sprintf('Transfer Orbit from %s to %s', strs{planet1}, strs{planet2}));
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
axis equal;
hold off;  % Optionally, you can call hold off here to ensure no further data is added to the current figure
end
