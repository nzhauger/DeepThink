clc;clear;close all;
%% --- Load MOOS Data ---
[x,y,lat,long,depth,heading,N] = MOOSalogScraper();

%% --- AUV Animation Parameters ---
% North 0 Degrees, CW = +Rotation

% IVER4 size in (x/y axes)
AUV_L = 2.5;
AUV_W = 0.2286; 

xpatch = [x(1) - AUV_L/2 x(1) + AUV_L/2 x(1) + AUV_L/2 x(1) - AUV_L/2];
ypatch = [y(1) - AUV_W/2 y(1) - AUV_W/2 y(1) + AUV_W/2 y(1) + AUV_W/2];

% buoy info
buoy_x = [72 71 69 68 69 71];
buoy_y = [-52 -50 -50 -52 -54 -54];
buoy_pos = [70, -52]; % buoy center

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
theta = 55; % has to do with where the sub is looking
d = 100: -0.2: 5; % distance to buoy. Shouldn't this be based on distance to final point?
%d_off = 100;

%phi = -30;
phi = heading(1);

dist_to_buoy = sqrt((buoy_pos(1) - x(1))^2 + (buoy_pos(2) - y(1))^2); %analogous to d_off
vec_to_buoy = buoy_pos - [x(1), y(1)];
heading_vec = [cosd(heading(1)), sind(heading(1))];
angle_diff = acosd(dot(vec_to_buoy, heading_vec) / (norm(vec_to_buoy)*norm(heading_vec)));

%off_calc = d_ff * tand(phi);
off_calc = dist_to_buoy * tand(-angle_diff); 
offset_arr = [off_calc, 0, 0, 0];
ar_pos = [-0.5, 0.5] + [off_calc, off_calc]; %position of ar tag

%ar tag image
artag = imread('artag.png');
artag2 = image(ax2, 'CData', flipud(artag), 'XData', ar_pos, 'YData', [1, 2]);

% Buoy rectangles
rect1 = rectangle(ax2, 'Position', [-0.75, 0, 1.5, 3] + offset_arr);
rect2 = rectangle(ax2, 'Position', [-2, -1, 4, 1] + offset_arr);
set(rect1, 'Visible','off');
set(rect2, 'Visible','off');
set(artag2,   'Visible','off');

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


% 9 smooth variances + square-style sharp edges
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
    %% 1. Update AUV animation
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
%% 2. Update Buoy & AR tag -----
    dist_to_buoy = sqrt((buoy_pos(1) - x(k))^2 + (buoy_pos(2) - y(k))^2);
    vec_to_buoy = buoy_pos - [x(k), y(k)];
    heading_vec = [cosd(heading(k)), sind(heading(k))];
    % I don't think angle diff is working as intended. 
    % Is the equation correct?
    angle_diff = acosd(dot(vec_to_buoy, heading_vec) / (norm(vec_to_buoy)*norm(heading_vec)));
    show_camera = (dist_to_buoy < 20) && (abs(angle_diff) < 60) && (depth(k) < 5);
    if show_camera
        cam_offset = dist_to_buoy * tand(-angle_diff);
        vw = 2 * dist_to_buoy * sind(55);
        xmax = vw / 2; xmin = -xmax;
        ymax = (1.85)*xmax+2; ymin = -ymax+4;
        axis(ax2, [xmin, xmax, ymin, ymax]);
        set(rect1, 'Position', [-0.75,0,1.5,3] + [cam_offset 0 0 0], 'Visible','on');
        set(rect2, 'Position', [-2,-1,4,1] + [cam_offset 0 0 0], 'Visible','on');
        set(artag2,   'XData', [-0.5 0.5] + cam_offset, 'Visible','on');
    else
        set(rect1, 'Visible','off');
        set(rect2, 'Visible','off');
        set(artag2,   'Visible','off');
    end


%% 3. Update Depth View animation --------
    
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


%% Finish loop
    phi = heading(k);
    drawnow;
    %pause(0.01);

end
%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unicycle Model kicks in
dt = 0.01;
% t = 0:dt:200;
% N = length(t);
% stride = 10;
next_wp = buoy_pos;
wp_i = 1;
%making empty vectors
heading_uni = zeros(2,1); x_uni = heading_uni; y_uni = heading_uni;
x_d = heading_uni; y_d = heading_uni; x_e = heading_uni; y_e = heading_uni;
theta_e = heading_uni; x_dot = heading_uni; y_dot = heading_uni;
theta_dot = heading_uni; u = heading_uni;
%populating index 1 with last known values
heading_uni(1) = 90 - heading(end - 1); x_uni(1) = x(end - 1); y_uni(1) = y(end - 1);
x_d(1) = next_wp(1); y_d(1) = next_wp(2);
V = 1.0; b = 1; Kp = 1; % parameters

%% Animation Loop for Unicycle Model
while dist_to_buoy > 2.5
    %% 1. Update AUV Animation 
    xpatch = [-AUV_L/2 AUV_L/2 AUV_L/2 -AUV_L/2];
    ypatch = [-AUV_W/2 -AUV_W/2 AUV_W/2 AUV_W/2];
   
    xrot = xpatch;
    yrot = ypatch;

    thetaAUV = heading_uni(1);
    M1 = [cosd(thetaAUV) -sind(thetaAUV)];
    M2 = [sind(thetaAUV) cosd(thetaAUV)];
    for j = 1:4
        rot = [M1; M2] * [xpatch(j); ypatch(j)];
        xrot(j) = rot(1);
        yrot(j) = rot(2);
    end
    
    set(AUVpatch, 'XData', xrot + x_uni(1), 'YData', yrot + y_uni(1))

    %% 2. Update Buoy & AR Tag View
    axes(ax2)
    % --- Unicycle Model Movement ---
    x_dot(2) = V * cosd(heading_uni(1));
    y_dot(2) = V * sind(heading_uni(1));
    theta_dot(2) = b * u(1);
    heading_uni(2) = heading_uni(1) + theta_dot(2)*dt;
    x_uni(2) = x_uni(1) + x_dot(2)*dt;
    y_uni(2) = y_uni(1) + y_dot(2)*dt;

    % --- Waypoint Navigation ---
    dx = next_wp(1) - x_uni(2);
    dy = next_wp(2) - y_uni(2);
    % if sqrt(dx^2 + dy^2) < 10 && wp_i < size(next_wp, 1)
    %     wp_i = wp_i + 1;
    % end
    x_d(2) = next_wp(1); y_d(2) = next_wp(2);
    x_e(2) = x_d(2) - x_uni(2); y_e(2) = y_d(2) - y_uni(2);
    heading_des = atan2d(y_e(2), x_e(2));
    if heading_des - heading_uni(2) > 180, heading_des = heading_des - 360;
    elseif heading_des - heading_uni(2) < -180, heading_des = heading_des + 360; end
    theta_e(2) = heading_des - heading_uni(2);
    u(2) = Kp * theta_e(2);


    dist_to_buoy = sqrt((buoy_pos(1) - x_uni(2))^2 + (buoy_pos(2) - y_uni(2))^2);
    % if dist_to_buoy < 0 % Not needed with the while loop
    %     V = 0; break;
    % end

    vec_to_buoy = buoy_pos - [x_uni(2), y_uni(2)];
    heading_vec = [90 - cosd(heading_uni(2)), 90 - sind(heading_uni(2))];
    angle_diff = acosd(dot(vec_to_buoy, heading_vec) / (norm(vec_to_buoy)*norm(heading_vec)));
    %show_camera = (dist_to_buoy < 20) && (abs(angle_diff) < 60);
    show_camera = (dist_to_buoy < 20) && (abs(theta_e(2)) < 60);
    if show_camera
        %cam_offset = dist_to_buoy * tand(-angle_diff);
        cam_offset = dist_to_buoy * tand(-theta_e(2));
        vw = 2 * dist_to_buoy * sind(55);
        xmax = vw / 2; xmin = -xmax;
        ymax = (1.85)*xmax+2; ymin = -ymax+4;
        axis([xmin, xmax, ymin, ymax]);
        set(rect1, 'Position', [-0.75,0,1.5,3] + [cam_offset 0 0 0], 'Visible','on');
        set(rect2, 'Position', [-2,-1,4,1] + [cam_offset 0 0 0], 'Visible','on');
        set(artag2,   'XData', [-0.5 0.5] + cam_offset, 'Visible','on');
    else
        set(rect1, 'Visible','off');
        set(rect2, 'Visible','off');
        set(artag2,   'Visible','off');
    end

    %update vectors for next cycle
    x_uni(1) = x_uni(2); y_uni(1) = y_uni(2); 
    heading_uni(1) = heading_uni(2);
    u(1) = u(2);

end

