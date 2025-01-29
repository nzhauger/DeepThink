% Tan Nguyen 28 Jan 2025
% Nick Hauger 28 Jan 2025
% TODO: check each str (contains?) and if it doesn't have the right 
% axis parameter, remove from the vector

clear variables; clc; close all;

alog = readlines('newfile.alog');
alog = splitlines(alog);

N = length(alog);

time = zeros(N,1);
xline = zeros(N, 1); i=1;
yline = zeros(N, 1); j=1;
zline = zeros(N, 1); z=1;


for k = 5:N
    alog(k) = erase(alog(k),'uSimMarine');
end

for k = 1:N
    if contains(alog(k),'NAV_X')
        xline(i) = alog(k);
        i=i+1;
    elseif contains(alog(k),'NAV_Y')
        yline(j) = alog(k);
        j = j+1;
    elseif contains(alog(k),'NAV_Z')
        zline(z) = alog(k);
        z = z+1;
    end

    % produces a NaN for each non NAV_K line
    % NAV_X(k) = extractAfter(alog(k), "X                uSimMarine      ");
    % NAV_Y(k) = extractAfter(alog(k), "Y                uSimMarine      ");
    % NAV_Z(k) = extractAfter(alog(k), "Z                uSimMarine      ");
    
    % did not work, elements did not match on both sides
    % NAV_X(k) = extractBetween(alog(k), 'NAV_X                uSimMarine      ', ' ');
    % NAV_Y(k) = extractBetween(alog(k), 'NAV_Y                uSimMarine      ', ' ');
    % NAV_Z(k) = extractBetween(alog(k), 'NAV_Z                uSimMarine      ', ' ');

end

