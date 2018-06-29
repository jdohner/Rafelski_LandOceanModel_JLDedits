% CLM_preprocess.m
%
% june 28, 2018
% 
% julia dohner
%
% file to process data from CLM 4.5 output received from Danica Lombardozzi

%ncdisp('TRENDY2017_S3_LAND_USE_FLUX.nc')

% days since 1860-01-01 00:00:00
time = ncread('TRENDY2017_S3_LAND_USE_FLUX.nc','time');

% history time interval endpoints (?)
time_bounds = ncread('TRENDY2017_S3_LAND_USE_FLUX.nc','time_bounds');

% coordinate longitude (degrees east)
% fill value/missing value = 9.999999616903162e+35
% increases from 0 to 358.75 in increments of 1.25 degres
lon = ncread('TRENDY2017_S3_LAND_USE_FLUX.nc','lon');

% coordinate latitude (degrees north)
% increments of 0.9424 degres
% from -90 to 90
% fill value/missing value = 9.999999616903162e+35
lat = ncread('TRENDY2017_S3_LAND_USE_FLUX.nc','lat');

% Size: 288x192x1884
% Dimensions: lon,lat,time
% total C emitted from land cover conversion 
%(smoothed over the year) and wood and grain product pools 
%(NOTE: not a net value)
%
% gC/m^2/s
% fill/missing value = 9.999999616903162e+35
LAND_USE_FLUX = ncread('TRENDY2017_S3_LAND_USE_FLUX.nc','LAND_USE_FLUX');
