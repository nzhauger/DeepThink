clc;clear;close all;
%% --- Load MOOS Data ---
[x,y,lat,long,depth,heading,N] = MOOSalogScraper();
% lat(length(lat)) = [];
% long(length(long)) = [];

% --- AUV Animation Parameters ---
% North 0 Degrees, CW = +Rotation

% IVER4 size in (x/y axes)
AUV_L = 2.5;
AUV_W = 0.2286; 

xpatch = [x(1) - AUV_L/2 x(1) + AUV_L/2 x(1) + AUV_L/2 x(1) - AUV_L/2];
ypatch = [y(1) - AUV_W/2 y(1) - AUV_W/2 y(1) + AUV_W/2 y(1) + AUV_W/2];

buoy_x = [72 71 69 68 69 71];
buoy_y = [-52 -50 -50 -52 -54 -54];

% --- Create Main Figure ---
figure('Position', [100, 100, 1600, 800]);
%% --- Subplot 1: AUV Position ---
ax1 = subplot(1,3,1);
s = hgtransform;
AUVpatch = patch(ax1,'XData', xpatch, 'YData', ypatch, 'Parent', s, 'FaceColor', 'b');
buoypatch = patch(ax1, 'XData', buoy_x, 'YData', buoy_y, 'FaceColor', 'r', 'EdgeColor', 'k');

axis equal;
axis([(min(x) - 10) (max(x) + 10) (min(y) - 10) (max(y) + 10)]);

title('AUV Movement');
xlabel('X');
ylabel('Y');

%% --- Subplot 2: Buoy and AR Tag View ---
ax2 = subplot(1,3,2);
hold(ax2, 'on');
theta = 55; %has to do with where the sub is looking
d = 100: -0.2: 5;
d_off = 100;
%phi = -30;
phi = heading(1);

off_calc = d_off * tand(phi); %vector of ???
off_arr = [off_calc, 0, 0, 0];
ar_pos = [-0.5, 0.5] + [off_calc, off_calc]; %position of ar tag

%ar tag image
artag = imread('artag.png');
artag2 = image(ax2, 'CData', flipud(artag), 'XData', ar_pos, 'YData', [1, 2]);

% Buoy rectangles
rect1 = rectangle(ax2, 'Position', [-0.75, 0, 1.5, 3] + off_arr);
rect2 = rectangle(ax2, 'Position', [-2, -1, 4, 1] + off_arr);

title('Buoy and AR Tag View');
xlabel('X');
ylabel('Y');

%% --- Subplot 3: Depth View ---
% Load depth vector and clean it
depth_vector = depth;

% Settings
n = length(depth_vector);
AUV_center_x = 50;                  % Fixed horizontal center
AUV_length = 4;                     % AUV body length
AUV_height = 1;                     % AUV body height

% Depth range
min_depth = 0;
max_depth = 30;
bottom_y = max_depth;

% Create figure
ax3 = subplot(1,3,3);
hold(ax3,'on');
axis equal;
xlim([AUV_center_x-10, AUV_center_x+10]);
ylim([min_depth - 2, max_depth]);
xlabel('Horizontal Distance (m)');
ylabel('Depth (m)');
title('AUV Side View - Depth Animation');
set(ax3, 'YDir', 'reverse');  % Depth increases downward

% Seabed (gravel/sand/rocky) with variation
gravel_x = linspace(AUV_center_x-10, AUV_center_x+10, 200);


%% 9 smooth variances + square-style sharp edges
smooth_noise = ...
    0.1 * sin(0.2 * gravel_x) + ...
    0.05 * cos(0.5 * gravel_x + pi/4) + ...
    0.07 * sin(1.0 * gravel_x + rand) + ...
    0.03 * cos(0.8 * gravel_x + rand) + ...
    0.06 * sin(0.4 * gravel_x + pi/2) + ...
    0.04 * cos(0.9 * gravel_x + pi/3) + ...
    0.02 * sin(1.3 * gravel_x + pi/6) + ...
    0.05 * sin(0.7 * gravel_x + 0.5) + ...
    0.03 * cos(1.1 * gravel_x + 0.7);

square_noise = ...
    0.1 * sign(sin(0.4 * gravel_x)) + ...
    0.05 * sign(sin(0.9 * gravel_x + pi/3));

rand_variation = 0.03 * randn(1, length(gravel_x));
gravel_y = bottom_y - 0.5 + smooth_noise + square_noise + rand_variation;

seabed = fill([gravel_x fliplr(gravel_x)], ...
              [gravel_y bottom_y*ones(size(gravel_x))], ...
              [0.76 0.69 0.57], 'EdgeColor', 'none'); % Sand color

% Water surface waves
wave_x = linspace(AUV_center_x-10, AUV_center_x+10, 100);
wave_line = plot(wave_x, 0.3*sin(0.3*wave_x), 'b', 'LineWidth', 2);

% AUV rectangle shape
depth_patch = rectangle('Position', ...
    [AUV_center_x - AUV_length/2, depth_vector(1) - AUV_height/2, AUV_length, AUV_height], ...
    'FaceColor', 'y', 'EdgeColor', 'k', 'Curvature', [0.3 0.3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Animation Loop (Synchronizing Both) ---

for k = 1:N-2
    %axes(ax1); %you dont need this if you specify patch is inside ax1
    %% Update AUV animation
    xpatch = [-AUV_L/2 AUV_L/2 AUV_L/2 -AUV_L/2];
    ypatch = [-AUV_W/2 -AUV_W/2 AUV_W/2 AUV_W/2];
   
    xrot = xpatch;
    yrot = ypatch;

    thetaAUV = 90 - heading(k);
    M1 = [cosd(thetaAUV) -sind(thetaAUV)];
    M2 = [sind(thetaAUV) cosd(thetaAUV)];
    for j = 1:4
        rot = [M1; M2] * [xpatch(j); ypatch(j)];
        xrot(j) = rot(1);
        yrot(j) = rot(2);
    end
    
    set(AUVpatch, 'XData', xrot + x(k), 'YData', yrot + y(k))
    
%% Update Depth View animation --------
    
    % Animation loop
    axes(ax3)
    % Update wave animation
    wave_line.YData = 0.3*sin(0.3*wave_x - 0.1*k);
    
    % Wiggle seabed with updated shape
    gravel_x_shifted = gravel_x + 0.1 * k;
    
    smooth_noise = ...
        0.1 * sin(0.2 * gravel_x_shifted) + ...
        0.05 * cos(0.5 * gravel_x_shifted + pi/4) + ...
        0.07 * sin(1.0 * gravel_x_shifted + rand) + ...
        0.03 * cos(0.8 * gravel_x_shifted + rand) + ...
        0.06 * sin(0.4 * gravel_x_shifted + pi/2) + ...
        0.04 * cos(0.9 * gravel_x_shifted + pi/3) + ...
        0.02 * sin(1.3 * gravel_x_shifted + pi/6) + ...
        0.05 * sin(0.7 * gravel_x_shifted + 0.5) + ...
        0.03 * cos(1.1 * gravel_x_shifted + 0.7);
    
    square_noise = ...
        0.1 * sign(sin(0.4 * gravel_x_shifted)) + ...
        0.05 * sign(sin(0.9 * gravel_x_shifted + pi/3));
    
    rand_variation = 0.03 * randn(1, length(gravel_x));
    gravel_y = bottom_y - 0.5 + smooth_noise + square_noise + rand_variation;
    seabed.YData = [gravel_y bottom_y*ones(size(gravel_x))];
    
    % Update AUV vertical position
    current_depth = depth_vector(k);
    depth_patch.Position(2) = current_depth - AUV_height/2;
    rotate(depth_patch,[1,0,0],30);

%% Update Buoy & AR Tag View
    axes(ax2)

    vw = 2*d(k)*sind(theta);
    xmax = vw/2;
    xmin = -xmax;
    ymax = (3/4)*xmax + 2;
    ymin = -ymax + 2;
    axis(ax2, [xmin, xmax, ymin, ymax]);

    %off_calc = d_off * tand(phi);
    off_calc = d_off * tand(heading(k));
    off_arr = [off_calc, 0, 0, 0];
    ar_pos = [-0.5, 0.5] + [off_calc, off_calc];

    if isvalid(rect1)
        set(rect1, "Position", [-0.75, 0, 1.5, 3] + off_arr);
    end
    if isvalid(rect2)
        set(rect2, "Position", [-2, -1, 4, 1] + off_arr);
    end
    if isvalid(artag2)
        set(artag2, 'XData', ar_pos);
    end

    % Update phi
    % Kp = 0.016;
    % phi = phi - Kp * phi;
    % phi = max(min(phi, 55), -55);

    drawnow;
    %pause(0.01);

end

