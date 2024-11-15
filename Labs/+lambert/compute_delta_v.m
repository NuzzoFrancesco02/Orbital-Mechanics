function d_v_lambert = compute_delta_v(t1, t2, planet1, planet2, mu)
% compute_delta_v - Computes the delta-v (change in velocity) required for a Lambert transfer
% between two planets in the solar system based on departure and arrival times.
%
% PROTOTYPE:
% d_v_lambert = compute_delta_v(t1, t2, planet1, planet2, mu)
%
% INPUT:
% t1 [1]          - Departure time in Modified Julian Date (MJD2000) [T]
% t2 [1]          - Arrival time in Modified Julian Date (MJD2000) [T]
% planet1 [int]   - Departure planet identifier (1 to 9 corresponding to planets) [-]
% planet2 [int]   - Arrival planet identifier (1 to 9 corresponding to planets) [-]
% mu [1]          - Gravitational parameter [L^3/T^2]
%
% OUTPUT:
% d_v_lambert [1] - The total delta-v required for the transfer [L/T]
%                   (sum of the required delta-vs at both departure and arrival)
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------

try


    % Retrieve the Keplerian elements and velocities of the departure and arrival planets
    kep1 = uplanet(t1, planet1);  % Keplerian elements of the departure planet at time t1
    kep2 = uplanet(t2, planet2);  % Keplerian elements of the arrival planet at time t2
    [r1, v1] = coord.kep2car_theta(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), mu);  % Position and velocity of departure planet
    [r2, v2] = coord.kep2car_theta(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6), mu);  % Position and velocity of arrival planet

    % Time of Flight (ToF) in seconds
    ToF = (t2 - t1) * 24 * 3600;  % Convert time difference from days to seconds

    % Lambert's problem solver with default settings
    Nrev = 0;              % Number of revolutions
    optionsLMR = 2;        % Lambert's method option (type of orbit)
    orbitType = 0;         % Orbit type (1: Prograde, -1: Retrograde)
    Ncase = 0;             % Lambert's case number (0: two-body transfer)
    
    % Solve Lambert's problem to find the velocity vectors at departure and arrival
    [at, pt, et, ~, vt1, vt2] = lambert.lambertMR(r1, r2, ToF, mu, orbitType, Nrev, Ncase, optionsLMR);
    
    % Calculate the delta-v as the sum of the velocity differences at departure and arrival
    d_v_lambert = norm(vt1 - v1) + norm(v2 - vt2);  % Total delta-v

catch
    % If an error occurs (e.g., invalid input or orbit calculation failure), return NaN
    d_v_lambert = NaN;
    disp('Error computing delta-v. Ensure correct input and orbital data.');
end

end
