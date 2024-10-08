function Earth3d(Rt, fg)
% Earth3d: Creates a 3D plot of Earth with a texture-mapped surface.
%
% This function plots a 3D model of Earth, optionally allowing you to set 
% the Earth's radius and figure background color.
%
% PROTOTYPE:
% Earth3d(Rt, fg)
%
% INPUT:
% Rt [1]        - Earth's radius in kilometers (default is 6378 km) [km]
% fg [1x1]      - Figure handle (optional) to set the figure background color
%
% OUTPUT:
% A 3D plot of Earth with a texture-mapped surface and customizable background.
%
% CONTRIBUTORS:
% Francesco Nuzzo
%
% VERSION HISTORY:
% 2024-10-8: First version
%
% -------------------------------------------------------------------------

    % Set default radius if not provided
    if nargin == 0
        Rt = 6378;  % Default Earth's radius in kilometers
    end
    
    % Specify the image used for the texture map of Earth
    Earth_image = "8081_earthmap10k.jpg";  % High-resolution Earth texture

    %% Figure Setup
    
    % Choose the color of the figure background
    background_plot = 'k';  % Set the background to black
    
    % Create the figure or set the background color if a figure handle is provided
    fg.Color = background_plot;  % Set the figure background color
    hold("on");  % Hold the current plot for adding more elements
    
    % Set the axes to equal scale to avoid distortion of Earth's shape
    axis("equal");
    
    % Label the axes for reference
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    
    % Set the initial view angle for the 3D plot (azimuth, elevation)
    view(120, 30);

    %% Create Earth Surface as a Wireframe Sphere
    
    % Define the number of panels for the sphere's surface resolution
    npanels = 100;  % More panels = smoother sphere surface
    
    % Create a 3D mesh grid for the ellipsoid (Earth) centered at (0,0,0)
    [x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, npanels);
    
    % Create a 3D plot of the Earth as a surface (initially wireframe, no color)
    globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');
    
    %% Apply Texture Map to the Globe
    
    % Load the image used as texture map for Earth's surface
    cdata = imread(Earth_image);
    
    % Set transparency level for the globe (1 = opaque, 0 = fully transparent)
    alpha = 1;
    
    % Apply the texture map to the surface using 'FaceColor' set to 'texturemap'
    % The 'CData' property is used to apply the image as the surface texture
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    
    %% Set the Grid and Background Colors
    
    % Set the axes grid color to white and increase grid line width
    ax = gca;  % Get the current axes handle
    set(ax, 'GridColor', [1, 1, 1], 'GridLineWidth', 2);  % White grid lines
    
    % Set the axes background color to black for better contrast
    set(ax, 'Color', 'k');
end
