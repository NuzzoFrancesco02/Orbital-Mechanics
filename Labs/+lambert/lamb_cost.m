function [d_v_min, t_dep, t_arr] = lamb_cost(t1_start, t1_end, t2_start, t2_end, planet1, planet2, mu, v_lim)
% lamb_cost - Computes the optimal transfer between two planets using Lambert's problem.
%
% PROTOTYPE:
% [d_v_min, t_dep, t_arr] = lamb_cost(t1_start, t1_end, t2_start, t2_end, planet1, planet2, mu)
%
% INPUT:
% t1_start [str]   - Start of departure window (datetime format) [T]
% t1_end   [str]   - End of departure window (datetime format) [T]
% t2_start [str]   - Start of arrival window (datetime format) [T]
% t2_end   [str]   - End of arrival window (datetime format) [T]
% planet1  [int]   - Identifier of departure planet (1: Mercury, ..., 9: Pluto) [-]
% planet2  [int]   - Identifier of arrival planet (1: Mercury, ..., 9: Pluto) [-]
% mu       [1]     - Gravitational parameter of the Sun [L^3/T^2]
%
% OUTPUT:
% d_v_min [1]      - Minimum Delta V required for transfer [L/T]
% t_dep   [str]    - Optimal departure date [T]
% t_arr   [str]    - Optimal arrival date [T]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
% 2024-11-15: Improved comments, plot clarity, and added date labels to axes.
%
% -------------------------------------------------------------------------

% Convert input dates to Modified Julian Date 2000 (MJD2000)
t1_start = timeConversion.date2mjd2000(t1_start);
t1_end = timeConversion.date2mjd2000(t1_end);
t2_start = timeConversion.date2mjd2000(t2_start);
t2_end = timeConversion.date2mjd2000(t2_end);

% Generate grids of departure and arrival dates
d_t1 = linspace(t1_start, t1_end, 300); % Departure times
d_t2 = linspace(t2_start, t2_end, 300); % Arrival times

% Initialize Lambert solver options
Nrev = 0;           % Number of revolutions
optionsLMR = 0;     % Algorithm option for Lambert solver
orbitType = 0;      % Prograde orbit
Ncase = 0;          % Single solution case

% Initialize Delta V matrix
d_v = NaN(length(d_t1), length(d_t2));

% Loop through all combinations of departure and arrival dates
for j = 1:length(d_t1)
    % Calculate departure planet position and velocity
    [kep1, ksun1] = uplanet(d_t1(j), planet1);
    [r1, v1] = coord.kep2car_theta(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), ksun1);

    for k = 1:length(d_t2)
        % Calculate arrival planet position and velocity
        [kep2, ~] = uplanet(d_t2(k), planet2);
        [r2, v2] = coord.kep2car_theta(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6), ksun1);

        % Time of flight (seconds)
        ToF = (d_t2(k) - d_t1(j)) * 24 * 60 * 60;

        % Solve Lambert's problem for the transfer
        [~, ~, ~, ~, vt1, vt2, ~, ~] = lambert.lambertMR(r1, r2, ToF, ksun1, orbitType, Nrev, Ncase, optionsLMR);

        % Calculate Delta V
        d_v1 = vt1 - v1; % Delta V at departure
        d_v2 = v2 - vt2; % Delta V at arrival
        d_v(j, k) = norm(d_v1) + norm(d_v2);

        % Discard high Delta V values
        if d_v(j, k) > v_lim
            d_v(j, k) = NaN;
        end
    end
end

% Find minimum Delta V and corresponding times
[d_v_dim1, pos1] = min(d_v, [], 1); % Minimum Delta V for each arrival date
[~, pos2] = min(d_v_dim1, [], 2);   % Overall minimum Delta V

t_dep = d_t1(pos1(pos2)); % Optimal departure date
t_start = d_t2(pos2);     % Optimal arrival date

% Fine-tune using optimization
objective = @(t) lambert.compute_delta_v(t(1), t(2), planet1, planet2, mu);
initial_guess = [t_dep, t_start];
options = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton', ...
                       'StepTolerance', 1e-16, 'OptimalityTolerance', 1e-16,'MaxIterations',1e8);

[optimal_times, d_v_min] = fminunc(objective, initial_guess, options);

% Convert optimized times to calendar dates
t_dep = timeConversion.mjd20002date(optimal_times(1));
t_arr = timeConversion.mjd20002date(optimal_times(2));

% Plot Delta V map
%figure
%surface(d_t1,d_t2,d_v')
figure;
levels = 5:1:10; % Define contour levels (only integers from 5 to 10)
contourf(d_t1, d_t2, d_v', levels, 'LineColor', 'none'); % Use filled contour for better visualization
hold on;
xlabel('Departure Date');
ylabel('Arrival Date');
title('Delta V Contour Map');
grid on;

% Configure colorbar
c = colorbar;
clim([min(min(d_v)) max(max(d_v))]); % Limit colorbar range
ylabel(c, '\DeltaV [km/s]', 'FontSize', 10);

% Convert MJD2000 to readable dates for axis labels
xticks = linspace(d_t1(1), d_t1(end), 6); % Select 6 points for x-axis labels
xticklabels = arrayfun(@(t) datestr(timeConversion.mjd20002date(t), 'dd-mmm-yyyy'), xticks, 'UniformOutput', false);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45);

yticks = linspace(d_t2(1), d_t2(end), 6); % Select 6 points for y-axis labels
yticklabels = arrayfun(@(t) datestr(timeConversion.mjd20002date(t), 'dd-mmm-yyyy'), yticks, 'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ytickangle(45);

% Highlight the optimal solution
plot(optimal_times(1), optimal_times(2), 'ro', 'MarkerFaceColor', 'r');
text(optimal_times(1), optimal_times(2), sprintf(' Min \\DeltaV = %.2f', d_v_min), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10,'Color','white');

% Display results
disp('Optimal Departure Date: ');
disp(datetime(t_dep));
disp('Optimal Arrival Date: ');
disp(datetime(t_arr));
disp(['Minimum Delta V: ', num2str(d_v_min)]);

end
