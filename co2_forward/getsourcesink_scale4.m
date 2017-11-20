% Get relevant land use and Joos model results
%
% Updates
%
% 4/13/07 - Lauren Elmegreen - use land use change 
% emissions that are extrapolated to the present by keeping values 
% constant at the last value, instead of increasing by 1.4%
% 10/25/12 - Lauren Rafelski - added annotations
% Oct 24, 2017 - JLD removed Aoc from output
%
% This file is for the record between 1800 and 2006+10/12

function [landusemo,ff1,fas,extratrop_landmo] = getsourcesink_scale4;

load joos_hilda_2011.mat
load fossil.mat
landnowppm_JLD = csvread('landnowppm_JLD.csv');
extratrop_landppm_JLD = csvread('extratrop_landppm_JLD.csv');


% interpolate land use changes to monthly resolution

%month = 1850:(1/12):2006; %old code
%month = 1800:(1/12):2010; % JLD on Oct 31, 2017 to extend land data to match time frame of ocean
month = 1800:(1/12):(2006+10/12);

%landnowppm begins as a 157x12 vector, gets interpolated to 
%landnowppm(:,1) contains the sample points
%landnowppm(:,12) contains the corresponding values
%month contains the coordinates of the query points
%month is a long vector (1x1873)
%landmonth = interp1(landnowppm(:,1),landnowppm(:,12),month);
landmonth = interp1(landnowppm_JLD(:,1),landnowppm_JLD(:,12),month);

%month2 = 1850:(1/12):2000; % old code
%month2 = 1800:(1/12):2010; % JLD on Oct 31, 2017
month2 = 1800:(1/12):(2006+10/12);
extralandmonth = interp1(extratrop_landppm_JLD(:,1),extratrop_landppm_JLD(:,2),month2);

landusemo(:,1) = month;
landusemo(:,2) = landmonth; % value in ppm

extratrop_landmo(:,1) = month2;
extratrop_landmo(:,2) = extralandmonth; % value in ppm

% other variables should just be loaded and passed