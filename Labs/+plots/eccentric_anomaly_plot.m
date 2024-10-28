function eccentric_anomaly_plot(a, e, mu, N, k, E0, t0)
% eccentric_anomaly_plot - Plots the eccentric anomaly as a function of normalized time for multiple eccentricities.
%
% PROTOTYPE:
%   eccentric_anomaly_plot(a, e, mu, N, k, E0, t0)
%
% DESCRIPTION:
%   This function calculates and plots the eccentric anomaly (E) over time for orbits with varying eccentricities.
%   The eccentric anomaly is solved for using Kepler's equation, and the results are displayed for different values
%   of eccentricity in a single plot.
%
% INPUT:
%   a       [1x1]    - Semi-major axis of the orbit                  [km]
%   e       [1xM]    - Array of eccentricities to plot               [-]
%   mu      [1x1]    - Gravitational parameter                       [km^3/s^2]
%   N       [1x1]    - Number of time steps for the plot             [-]
%   k       [1x1]    - Number of orbital periods to simulate         [-]
%   E0      [1x1]    - Initial eccentric anomaly                     [rad]
%   t0      [1x1]    - Initial time                                  [s]
%
% OUTPUT:
%   None. The function generates a plot of the eccentric anomaly over normalized time.
%
% CONTRIBUTORS:
%   Francesco Nuzzo
%
% VERSIONS:
%   2024-10-10: First version
%
% -------------------------------------------------------------------------

    % Calculate the orbital period
    T = 2 * pi * sqrt(a^3 / mu); % Orbital period [s]
    
    % Generate a time vector from the initial time to the end of k orbital periods
    t = linspace(t0, k * T, N);  % Time vector [s]

    % Initialize a legend array for plot labeling
    legend_names = strings(1, length(e));
    
    % Loop over each eccentricity in the input array
    for i = 1:length(e)
        % Solve Kepler's equation for the eccentric anomaly (E) over time
        [~, E] = elliptic.kepler_equation(t, e(i), a, mu, t0, E0);
        
        % Plot the eccentric anomaly in degrees, normalized by orbital period (t/T)
        plot(t ./ T, rad2deg(E), 'LineWidth', 2); hold on; grid on;
        
        % Create a legend entry for each eccentricity value
        legend_names(i) = ['e = ' num2str(e(i))];
    end
    
    % Add legend, set axis limits and labels, and format plot
    legend(legend_names);  % Display eccentricities in legend
    xlim([0, t(end) / T]); % Limit x-axis to full orbital period range
    ylim([0, max(rad2deg(E))]); % Set y-axis limit based on max eccentric anomaly
    yticks(linspace(0, max(rad2deg(E)), 7)); % Define y-axis tick spacing
    
    % Set plot title and axis labels
    title('Eccentric Anomaly vs. Time'); 
    xlabel('t [T]'); ylabel('E [deg]');
end
