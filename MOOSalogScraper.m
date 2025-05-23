function [lat,long,depth,heading,N] = MOOSalogScraper()
% Tan Nguyen 28 Jan 2025
% Nick Hauger 23 MAY 2025
% this one works

clear variables; clc; close all;

alog = readlines('nav_depth.alog'); %this is just to grab N rows
alog = splitlines(alog);

%N = length(alog);
% alog = alog(1:N);
N = length(alog);

time = zeros(N,1);
lat = zeros(N,1);
long = zeros(N,1);
depth = zeros(N,1);
heading = zeros(N,1);

%creating vectors of the individual log files
latlog = readlines('nav_lat.alog');
latlog = splitlines(latlog);
longlog = readlines('nav_long.alog');
longlog = splitlines(longlog);
depthlog = readlines('nav_depth.alog');
depthlog = splitlines(depthlog);
hdglog = readlines('nav_heading.alog');
hdglog = splitlines(hdglog);

%grabbing values
for k = 1:N-1
    time(k) = extractBefore(latlog(k),',NAV_LAT');
    lat(k) = extractAfter(latlog(k),'NAV_LAT,');
    long(k) = extractAfter(longlog(k),'NAV_LONG,');
    depth(k) = extractAfter(depthlog(k),'NAV_DEPTH,');
    heading(k) = extractAfter(hdglog(k),'NAV_HEADING,');
end

%lat_deg_per_meter = 1 / 111000;             % Approximate
%lon_deg_per_meter = 1 / (85000);            % At ~44°N
%lat = lat./111000;
%long = long ./(cosd(43.825)*111000);

