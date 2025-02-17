% Tan Nguyen 28 Jan 2025
% Nick Hauger 28 Jan 2025

clear variables; clc; close all;

alog = readlines('newfile.alog');
alog = splitlines(alog);

N = length(alog);
alog = alog(5:N);
N = length(alog);

time = zeros(N,1);
xval = zeros(N,1);
yval = zeros(N,1);
zval = zeros(N,1);

for k = 1:N
    alog(k) = erase(alog(k),'uSimMarine');
    alog(k) = erase(alog(k),' ');
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

%trying to plot this 
figure(1);
hold on;
axis([-150 100 -150 100])
AUV_trail = plot3([xval(1) xval(2)],[yval(1) yval(2)],[zval(1) zval(2)],'--'); 
shipbody_x = [-1, 0 1 -1];
shipbody_y = [0, -1 0 0];
shipbody_z = [0, 0, 0, 0];
s = hgtransform;
patch('XData', shipbody_x, 'YData', shipbody_y, 'ZData', shipbody_z, 'Parent',s)

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
