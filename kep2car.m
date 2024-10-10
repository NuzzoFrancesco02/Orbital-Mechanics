function [r,v,theta] = kep2car(a, e, i, OM, om, theta, mu)
%   [r,v] = kep2car(a,e,i,OM,w,theta,mu)
%
% DESCRIPTION
%   Conversion from Keplerian elements to Cartesian coordinates. It works
%   also with vectorial theta [1xN]
%
% INPUT
%   a      [1x1]   Semi major axis                         [km]
%   e      [1x1]   Eccentricity                            [-]
%   i      [1x1]   Inclination                             [rad]
%   OM     [1x1]   Right Ascension of the Ascending Node   [rad]
%   w      [1x1]   Pericentre anomaly                      [rad]
%   theta  [1xN]   True anomaly                            [rad]
%   mu     [1x1]   Gravitational parameter                 [km^3/s^2]
%
% OUTPUT
%   r   [3xN]   Position vector             [km]
%   v   [3xN]   Velocity vector             [km/s]
%   theta [1xN] True anomaly (useful when e >= 1) [rad]

p = a*(1-e^2);
if e == 1
    p = 2*a;
    rp_max = 1e5; 
    if length(theta) ~= 1
        theta = -floor(acos(p/rp_max-1)):0.005:floor(acos(p/rp_max-1));
    end
elseif e > 1
    rp_max = 1e5; 
    if length(theta) ~= 1
        theta = -floor(acos((p/rp_max-1)/e)):0.01:floor(acos((p/rp_max-1)/e));
    end
end
l = length(theta);
R3 = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R1 = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
R = R1*R2*R3;
R = R';
r = []; v = [];
j = 0; 

while j < l 
    j = j + 1;
    if e < 1
        r_pf = ((p./(1+e.*cos(theta(j)))).*[cos(theta(j)) sin(theta(j)) 0]');
        v_pf = (sqrt(mu/p).*[-sin(theta(j)) e+cos(theta(j)) 0]');
    elseif e > 1
        r_pf = real((p./(1+e.*cos(theta(j)))).*[cos(theta(j)) sin(theta(j)) 0]');
        v_pf = imag(sqrt(mu/p).*[-sin(theta(j)) e+cos(theta(j)) 0]');
    elseif e == 1
         r_pf = ((p./(1+e.*cos(theta(j)))).*[cos(theta(j)) sin(theta(j)) 0]');
         v_pf = (sqrt(mu/p).*[-sin(theta(j)) e+cos(theta(j)) 0]');
    end
    r = [r; (R*r_pf)'];
    v = [v; (R*v_pf)'];
end
 