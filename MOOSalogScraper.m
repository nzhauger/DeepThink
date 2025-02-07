% Tan Nguyen 28 Jan 2025
% Nick Hauger 28 Jan 2025

clear variables; clc; close all;

% read the alog file, split it at the newlines
alog = readlines('newfile.alog');
alog = splitlines(alog);

% get rid of the first 4 lines that don't have data.
N = length(alog);
alog = alog(5:N);
N = length(alog);

% preallocate arrays
time = zeros(N,1);
xval = zeros(N,1);
yval = zeros(N,1);
zval = zeros(N,1);

% This is the loop where we go through the alog lines and extract the data into its particular arrays. 
for k = 1:N
    alog(k) = erase(alog(k),'uSimMarine'); % erasing the uSimMarines in the line. We know its coming from there.
    alog(k) = erase(alog(k),' '); % remove spaces
    % check for the variables we're looking for in the line and grab the data.
    if contains(alog(k),'NAV_X')
        time(k) = extractBefore(alog(k),'NAV_X');
        xval(k:N) = extractAfter(alog(k),'NAV_X');
    elseif contains(alog(k),'NAV_Y')
        time(k) = extractBefore(alog(k),'NAV_Y');
        yval(k:N) = extractAfter(alog(k),'NAV_Y');
    elseif contains(alog(k),'NAV_Z')
        time(k) = extractBefore(alog(k),'NAV_Z');
        zval(k:N) = extractAfter(alog(k),'NAV_Z');
    end
end

% trying to plot this 
figure(1); 
hold on; 
axis([-1000 4000 -1000 3000]) 
AUV_trail = plot3([xval(1) xval(2)],[yval(1) yval(2)],[zval(1) zval(2)],'--'); % plotting the first point of a trail for the UUV
% creating polygon points for a ship body. Right now its a triangle (first and last points should complete a shape)
shipbody_x = [-1, 0 1 -1]; 
shipbody_y = [0, -1 0 0]; 
shipbody_z = [0, 0, 0, 0]; 
s = hgtransform;
patch('XData', shipbody_x, 'YData', shipbody_y, 'ZData', shipbody_z, 'Parent',s) % create the patch object to plot it

% plotting 'over time' by translating the matrix and moving the trail. It draws as far as your computer runs it,
% not exactly per second.
for i = 1:N-1
    %plot3(xval(i),yval(i),zval(i),'*')
    s.Matrix = makehgtform('translate', [xval(i+1) yval(i+1) zval(i+1)]);
    set(AUV_trail, 'XData',xval(1:i))
    set(AUV_trail, 'YData',yval(1:i))
    set(AUV_trail, 'ZData',zval(1:i))
    xlabel('X position (m)'); ylabel('Y position (m)'); zlabel('Depth (m)');
    %hold off;
    drawnow;
end
