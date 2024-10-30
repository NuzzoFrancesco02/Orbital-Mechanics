function static_orbit(r)
% This function plots the 3D trajectory of an orbit over time.
% The Earth is displayed in 3D, and the orbit is plotted as a line.
%
% INPUT:
% r [nx3]    - Position vector at each time step (x, y, z) [km]
%
% OUTPUT:
% A 3D plot of the orbit with the Earth centered at the origin.
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSIONS:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

    if ~exist("figure")
        figure()
    end
    
    % Call the function Earth3d to display a 3D model of the Earth
    plots.Earth3d 
    
    % Plot the 3D orbit path using position vectors
    [~,min_dim] =  min(size(r));
    if  min_dim == 1
        r = r';
    end
    plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 2); 
    
    % Label the axes for better understanding of the plot
    xlabel('X [km]'); 
    ylabel('Y [km]'); 
    zlabel('Z [km]');
    
    % Add a title to the plot to describe the content
    title('2BP-ode Orbit');
    
    % Ensure that the axes are equally scaled to avoid distortion
    axis equal;
    
    % Turn on the grid to make it easier to visually interpret the plot
    grid on;
end
