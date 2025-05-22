clear variables; clc; close all;

% Set figure window positions to left and right halves of the screen
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');  % [left bottom width height]
screen_width = screenSize(3);
screen_height = screenSize(4);

% Load waypoint and depth data
[x_a, y_a, z_a, alog_N] = MOOSalogScraper();
load('depthval.mat');
depthval2 = nonzeros(depthval);
depth_vector = depthval2(:);
depth_vector = depth_vector(isfinite(depth_vector));

% --- Parameters ---
dt = 0.01;
t = 0:dt:200;
N = length(t);
stride = 10;
next_wp = [x_a(1:stride:end), y_a(1:stride:end)];
wp_i = 1;
theta = zeros(N,1); x = theta; y = theta;
x_d = theta; y_d = theta; x_e = theta; y_e = theta;
theta_e = theta; x_dot = theta; y_dot = theta; theta_dot = theta; u = theta;
theta(1) = 270;
x(1) = x_a(1); y(1) = y_a(1);
x_d(1) = next_wp(1,1); y_d(1) = next_wp(1,2);
V = 20; b = 5; Kp = 1;

% --- Buoy Info ---
buoy_pos = [x_a(end), y_a(end)];
artag = imread('artag.png');

% === Figure 1: Top-Down View + Camera ===
f1 = figure(1); clf;
set(f1, 'Position', [50, 100, screen_width/2.6, screen_height/1.75]);  % Left half
subplot(1,2,1); hold on;
axis equal; xlim([0 50]); ylim([-100 10]);
ship_trail = plot([x(1) x(2)], [y(1) y(2)], '--');
current_waypoint = plot(x_d(1), y_d(1), 'x');
set(0, 'DefaultFigurePosition', [10, 50, 1500, 720])
s = hgtransform;
patch('XData', 0.5*[-4 4 4 -4]*1.3, 'YData', 0.5*[-0.6 -0.6 0.6 0.6]*1.3, 'Parent', s, 'FaceColor', 'b');
plot(buoy_pos(1), buoy_pos(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'LineWidth', 1.5);

subplot(1,2,2); hold on;
cam_rect1 = rectangle('Position',[-0.75,0,1.5,3], 'Visible','off'); 
cam_rect2 = rectangle('Position',[-2,-1,4,1], 'Visible','off'); 
cam_img   = image(flipud(artag), 'XData', [-0.5 0.5], 'YData', [1, 2], 'Visible','off');

% === Figure 2: Side-View Depth Animation ===
AUV_center_x = 50; AUV_length = 4; AUV_height = 1;
min_depth = 0; max_depth = 30; bottom_y = max_depth;
gravel_x = linspace(AUV_center_x-10, AUV_center_x+10, 200);
wave_x = linspace(AUV_center_x-10, AUV_center_x+10, 100);

f2 = figure(2); clf; hold on;
set(f2, 'Position', [screen_width/2, 100, screen_width/2.5, screen_height/1.5]);  % Right half
axis equal;
xlim([AUV_center_x-10, AUV_center_x+10]);
ylim([min_depth - 2, max_depth]);
xlabel('Horizontal Distance (m)');
ylabel('Depth (m)');
title('AUV Side View - Depth Animation');
set(gca, 'YDir', 'reverse');

gravel_x_shifted = gravel_x;
smooth_noise = 0.1 * sin(0.2 * gravel_x_shifted);  % Placeholder for structure
square_noise = 0.1 * sign(sin(0.4 * gravel_x_shifted));
rand_variation = 0.03 * randn(1, length(gravel_x));
gravel_y = bottom_y - 0.5 + smooth_noise + square_noise + rand_variation;

seabed = fill([gravel_x fliplr(gravel_x)], ...
              [gravel_y bottom_y*ones(size(gravel_x))], ...
              [0.76 0.69 0.57], 'EdgeColor', 'none');
wave_line = plot(wave_x, 0.3*sin(0.3*wave_x), 'b', 'LineWidth', 2);
AUV_patch = rectangle('Position', ...
    [AUV_center_x - AUV_length/2, depth_vector(1) - AUV_height/2, AUV_length, AUV_height], ...
    'FaceColor', 'y', 'EdgeColor', 'k', 'Curvature', [0.3 0.3]);

% === Unified Real-Time Simulation Loop ===
for i = 2:N
    % --- Unicycle Model Movement ---
    x_dot(i) = V * cosd(theta(i-1));
    y_dot(i) = V * sind(theta(i-1));
    theta_dot(i) = b * u(i-1);
    theta(i) = theta(i-1) + theta_dot(i)*dt;
    x(i) = x(i-1) + x_dot(i)*dt;
    y(i) = y(i-1) + y_dot(i)*dt;

    % --- Waypoint Navigation ---
    dx = next_wp(wp_i, 1) - x(i);
    dy = next_wp(wp_i, 2) - y(i);
    if sqrt(dx^2 + dy^2) < 10 && wp_i < size(next_wp, 1)
        wp_i = wp_i + 1;
    end
    x_d(i) = next_wp(wp_i, 1); y_d(i) = next_wp(wp_i, 2);
    x_e(i) = x_d(i) - x(i); y_e(i) = y_d(i) - y(i);
    theta_des = atan2d(y_e(i), x_e(i));
    if theta_des - theta(i) > 180, theta_des = theta_des - 360;
    elseif theta_des - theta(i) < -180, theta_des = theta_des + 360; end
    theta_e(i) = theta_des - theta(i);
    u(i) = Kp * theta_e(i);

    % --- Stop if Close to Buoy ---
    dist_to_buoy = sqrt((buoy_pos(1) - x(i))^2 + (buoy_pos(2) - y(i))^2);
    if dist_to_buoy < 5
        V = 0; break;
    end

    % === Update Top-Down and Camera View (Figure 1) ===
    figure(1);
    subplot(1,2,1);
    s.Matrix = makehgtform('translate', [x(i) y(i) 0], 'zrotate', deg2rad(theta(i)));
    set(ship_trail, 'XData', x(1:i), 'YData', y(1:i));
    set(current_waypoint, 'XData', x_d(i), 'YData', y_d(i));
    
    subplot(1,2,2);
    vec_to_buoy = buoy_pos - [x(i), y(i)];
    heading_vec = [cosd(theta(i)), sind(theta(i))];
    angle_diff = acosd(dot(vec_to_buoy, heading_vec) / (norm(vec_to_buoy)*norm(heading_vec)));
    show_camera = (dist_to_buoy < 20) && (abs(angle_diff) < 60);
    if show_camera
        cam_offset = dist_to_buoy * tand(-angle_diff);
        vw = 2 * dist_to_buoy * sind(55);
        xmax = vw / 2; xmin = -xmax;
        ymax = (1.85)*xmax+2; ymin = -ymax+4;
        axis([xmin, xmax, ymin, ymax]);
        set(cam_rect1, 'Position', [-0.75,0,1.5,3] + [cam_offset 0 0 0], 'Visible','on');
        set(cam_rect2, 'Position', [-2,-1,4,1] + [cam_offset 0 0 0], 'Visible','on');
        set(cam_img,   'XData', [-0.5 0.5] + cam_offset, 'Visible','on');
    else
        set(cam_rect1, 'Visible','off');
        set(cam_rect2, 'Visible','off');
        set(cam_img,   'Visible','off');
    end

    % === Update Side-View Depth (Figure 2) ===
    figure(2);
    wave_line.YData = 0.3 * sin(0.3 * wave_x - 0.1 * i);
    gravel_x_shifted = gravel_x + 0.1 * i;
    smooth_noise = 0.1 * sin(0.2 * gravel_x_shifted) + ...
                   0.05 * cos(0.5 * gravel_x_shifted + pi/4);  % Simplified
    square_noise = 0.1 * sign(sin(0.4 * gravel_x_shifted));
    rand_variation = 0.03 * randn(1, length(gravel_x));
    gravel_y = bottom_y - 0.5 + smooth_noise + square_noise + rand_variation;
    seabed.YData = [gravel_y bottom_y*ones(size(gravel_x))];
    
    if i <= length(depth_vector)
        current_depth = depth_vector(i);
        AUV_patch.Position(2) = current_depth - AUV_height/2;
    end

    pause(0.03);
    drawnow;
end
