function Earth3d(Rt,fg)
if nargin==0
    Rt=6378;
end
Earth_image = "8081_earthmap10k.jpg";
%Earth_image = "BlueMarble.png";
% Earth_image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

%% Figure

% Choose the color of the figure background
background_plot = 'k';

% Create the figure
fg.Color = background_plot;
%figure('Color', background_plot);
hold("on");
%grid("on");

% Set the axes scale equal
axis("equal");

% Put the axes labels
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

% Set initial view
view(120,30);

%% Create Earth surface as a wireframe

% Define the number of panels to be used to model the sphere 
npanels = 100;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, npanels);

% Create the globe with the surf function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

%% Texturemap the globe

% Load Earth image for texture map
cdata = imread(Earth_image);

% Set the transparency of the globe: 1 = opaque, 0 = invisible
alpha = 1; 

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded 
% from the image. Finally, set the transparency and remove the edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
ax = gca; set(ax,'Color','k');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.ZAxis.Visible = 'off';
end