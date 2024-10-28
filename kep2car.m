function [r,v,th] = kep2car(a, e, i, OM, om, th, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
%   a             [1x1]   Semi-major axis             [km]
%   e             [1x1]   Eccentricity                [-]
%   i             [1x1]   Inclination                 [rad]
%   OM            [1x1]   RAAN                        [rad]
%   om            [1x1]   Pericentre anomaly          [rad]
%   th            [1x1]   True anomaly                [rad]
%   mu            [1x1]   Gravitational parameter     [km^3/s^2]
%
% OUTPUT:
%   r             [3x1]   Position vector             [km]
%   v             [3x1]   Velocity vector             [km/s]

p = a*(1-e^2);
if e == 1
    p = 2*a;
    rp_max = 1e5; 
    th = -floor(acos(p/rp_max-1)):0.005:floor(acos(p/rp_max-1));
elseif e > 1
    rp_max = 1e5; 
    th = -floor(acos((p/rp_max-1)/e)):0.01:floor(acos((p/rp_max-1)/e));
end
l = length(th);
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
        r_pf = ((p./(1+e.*cos(th(j)))).*[cos(th(j)) sin(th(j)) 0]');
        v_pf = (sqrt(mu/p).*[-sin(th(j)) e+cos(th(j)) 0]');
    elseif e > 1
        r_pf = real((p./(1+e.*cos(th(j)))).*[cos(th(j)) sin(th(j)) 0]');
        v_pf = imag(sqrt(mu/p).*[-sin(th(j)) e+cos(th(j)) 0]');
    elseif e == 1
         r_pf = ((p./(1+e.*cos(th(j)))).*[cos(th(j)) sin(th(j)) 0]');
         v_pf = (sqrt(mu/p).*[-sin(th(j)) e+cos(th(j)) 0]');
    end
    %if min(r_pf)>-1e5 && max(r_pf)<1e5
        r = [r R*r_pf];
        v = [v R*v_pf];
    %end
    
end
 