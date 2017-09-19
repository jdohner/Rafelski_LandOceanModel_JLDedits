% file GetSourceSink.m
% formerly "getsourcesink_scale3.m"
%
% author Lauren Rafelski, modified by Julia Dohner
% 
% brief Get relevant land use and Joos model results
%
% changes by LR
% 4/13/07 - Lauren Elmegreen - use land use change emissions that are 
% extrapolated to the present by keeping values constant at the last value,
% instead of increasing by 1.4%
% 10/25/12 - Lauren Rafelski - added annotations
%
% TODO:
% finish briefs
% see if can rename  variables in joos_hilda_2011.m file to match my other 
% changes i.e. ff1 -> fossilFuelData, fas -> airSeaFlux

function [landUse,fossilFuelData,airSeaFlux,Aoc,extratropLandUse] = GetSourceSink;

% Get relevant land use data and Joos model results
load joos_hilda_2011.mat
load fossil.mat

% renaming outputs so match convention
fossilFuelData = ff1; % retrieved from joos_hilda_2011.mat
airSeaFlux = fas; % retrieved from joos_hilda_2011.mat

%% land use emissions

% create month array for all months between 1850 and 2006
landUseMonths = 1850:(1/12):2006; 
% interpolate land use emissions to monthly (in ppm/yr)
landUseData = interp1(landnowppm(:,1),landnowppm(:,12), landUseMonths);

% land use emissions output matrix
landUse(:,1) = landUseMonths;
landUse(:,2) = landUseData; % value in ppm

%% extra tropical land use emissions 

% create month array for all months between 1850 and 2000
extratropMonths = 1850:(1/12):2000;
% interpolate extratropical land use emissions to monthly (in ppm/yr)
extratropData = interp1(extratrop_landppm(:,1),extratrop_landppm(:,2),extratropMonths);

% extratropical land use emissions output matrix
extratropLandUse(:,1) = extratropMonths;
extratropLandUse(:,2) = extratropData; % value in ppm

