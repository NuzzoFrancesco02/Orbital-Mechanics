function [a, e, i, OM, om, th] = car2kep_theta(r, v, mu)
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
%   [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in % radians.
%
% INPUT:
%   r             [Nx3]   Position vector                 [km]
%   v             [Nx3]   Velocity vector                 [km/s]
%   mu            [1x1]   Gravitational parameter         [km^3/s^2]
%
% OUTPUT:
%   a             [1x1]   Semi-major axis                 [km]
%   e             [1x1]   Eccentricity                    [-]
%   i             [1x1]   Inclination                     [rad]
%   OM            [1x1]   RAAN                            [rad]
%   om            [1x1]   Pericentre anomaly              [rad]
%   th            [Nx1]   True anomaly                    [rad]

x = [1 0 0]';
y = [0 1 0]';
z = [0 0 1]';
l = max(size(r,1));


E = 0.5*norm(v)^2-mu/norm(r);
a = -mu/(2*E); %
e_vet = ((norm(v)^2-mu/norm(r))*r-dot(r,v)*v)/mu;
e = norm(e_vet); %
if e < 1e-15
    e = 0;
end

h = cross(r(1,:),v(1,:)); 
i = acos(dot(h,z)/norm(h));


N = cross(z,h);
if N == [0 0 0]
    N = x;
end

if e ~= 0
    if i ~= 0
        
        OM = acos(dot(N,x)/norm(N));
        if dot(N,y) <= 0
            OM = 2*pi - OM;
        end
    else
        OM = 0;
    end
   
    om = acos(dot(N,e_vet)/(norm(N)*e));
    if dot(e_vet,z) <= 0
        om = 2*pi - om;
    end

elseif e == 0 
    om = 0; 
    if i ~= 0
    
        OM = acos(dot(N,x)/norm(N));
        if dot(N,y) <= 0
            OM = 2*pi - OM;
        end
    else
        OM = 0;
    end
end
th = [];
for j = 1 : l
    vr = dot(r(j,:),v(j,:))./norm(r);
    if e == 0
        th = [th acos(dot(x,r(j,:))/norm(r(j,:)))];
    else
        if vr < 0
            th = [th 2*pi - real(acos(dot(e_vet,r(j,:))/(e*norm(r))))];
        else
            th = [th real(acos(dot(e_vet,r(j,:))/(e*norm(r))))];
        end
    end
end
