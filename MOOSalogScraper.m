% Tan Nguyen 28 Jan 2025
% Nick Hauger 28 Jan 2025
% TODO: check each str (contains?) and if it doesn't have the right 
% axis parameter, remove from the vector

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
for i = 1:N
    hold on;
    plot3(xval(i),yval(i),zval(i),'*')
    hold off;
    getframe;
end

