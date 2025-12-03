%% Load and visualize risk map
clear; clc; close all;

%% 1. Load data
risk = readmatrix('riskMap.csv');   % total fatality risk per cell
lat  = readmatrix('lat.csv');       % latitude of each cell
lon  = readmatrix('lon.csv');       % longitude of each cell

% Check sizes
if ~isequal(size(risk), size(lat), size(lon))
    error('riskMap, lat and lon must all have the same size!');
end

[nRows, nCols] = size(risk);
fprintf('Loaded risk map: %d rows x %d cols\n', nRows, nCols);

%% 2. Visualization using grid indices
figure;
imagesc(risk);
set(gca, 'YDir', 'normal');
axis equal tight;
colormap(jet);
colorbar;
title('Ground Risk Map (Grid Indices)');
xlabel('Column Index j');
ylabel('Row Index i');

%% 3. Visualization using geographic coordinates
figure;
pcolor(lon, lat, risk);
shading flat;
axis equal tight;
colormap(jet);
colorbar;
title('Ground Risk Map (Geo Coordinates)');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% 4. Compute grid resolution (meters per cell)
R = 6371000; % Earth radius

% Horizontal (1,1) -> (1,2)
dlat = deg2rad(lat(1,2) - lat(1,1));
dlon = deg2rad(lon(1,2) - lon(1,1));
latm = deg2rad((lat(1,1)+lat(1,2))/2);
dx = R * sqrt( (dlat)^2 + (cos(latm)*dlon)^2 );

% Vertical (1,1) -> (2,1)
dlat = deg2rad(lat(2,1) - lat(1,1));
dlon = deg2rad(lon(2,1) - lon(1,1));
latm = deg2rad((lat(1,1)+lat(2,1))/2);
dy = R * sqrt( (dlat)^2 + (cos(latm)*dlon)^2 );

fprintf('Grid resolution: dx = %.2f m, dy = %.2f m\n', dx, dy);
