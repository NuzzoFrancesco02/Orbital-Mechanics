function [alpha, delta, l, m, n] = eci2ra(r)
% eci2ra - Convert ECI coordinates to right ascension and declination.
%
% PROTOTYPE:
% [alpha, delta, l, m, n] = eci2ra(r)
%
% INPUT:
% r [nx3]    - Position vectors in the Earth-Centered Inertial (ECI) frame [L]
% 
% OUTPUT:
% alpha [nx1] - Right ascension (alpha) in radians [rad]
% delta [nx1] - Declination (delta) in radians [rad]
% l [nx1]     - Unit vector component in x-direction [unitless]
% m [nx1]     - Unit vector component in y-direction [unitless]
% n [nx1]     - Unit vector component in z-direction [unitless]
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-10: First version
%
% -------------------------------------------------------------------------
    
    % Check the dimensions of r and transpose if necessary
    if size(r, 2) > 3
        r = r'; % Transpose to ensure r is a column vector
    end
    
    % Calculate the norm (magnitude) of the position vector and unit components
    r_norm = vecnorm(r, 2, 2);  % Compute the Euclidean norm along rows
    l = r(:, 1) ./ r_norm; % Component in the x-direction (normalized)
    m = r(:, 2) ./ r_norm; % Component in the y-direction (normalized)
    n = r(:, 3) ./ r_norm; % Component in the z-direction (normalized)
    
    % Calculate delta (declination)
    delta = asin(n); % declination calculated from the z-component
    
    % Calculate alpha (right ascension)
    alpha = acos(l ./ cos(delta)); % Right ascension from the x-component
    
    % Correction for the quadrant
    alpha(m < 0) = 2 * pi - alpha(m < 0); % Adjust for negative y-components

    %alpha( abs( delta ) ~= pi/2 ) = 0;
    % Ensure alpha is in radians
    alpha = wrapTo2Pi(real(alpha)); % Force alpha to be in the range [0, 2*pi]
    
end
