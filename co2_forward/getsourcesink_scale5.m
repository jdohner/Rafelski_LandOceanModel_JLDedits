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

function [ff1,landusemo,extratrop_landmo] = getsourcesink_scale5(year);

load joos_hilda_2011.mat
landnowppm_JLD = csvread('landnowppm_JLD.csv');
extratrop_landppm_JLD = csvread('extratrop_landppm_JLD.csv');

% shortening ff vector to begin at start_year
ff1_start = find(ff1(:,1) == 1800);
ff1_end = find(ff1(:,1) == 2009+(7/12));
ff1 = ff1(ff1_start:ff1_end,:);

% interpolate land use changes to monthly resolution

month = 1800:(1/12):2006;
landmonth = interp1(landnowppm_JLD(:,1),landnowppm_JLD(:,12),month);

month2 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
extralandmonth = interp1(extratrop_landppm_JLD(:,1),extratrop_landppm_JLD(:,2),month2);

landusemo(:,1) = month;
landusemo(:,2) = landmonth; % value in ppm

extratrop_landmo(:,1) = month2;
extratrop_landmo(:,2) = extralandmonth; % value in ppm

last_ind_2006 = length(month);
last_ind_2000 = length(month2);

%Extend extratropical emissions by assuming emissions are zero
%extratrop_landmo(:,1) = landusemo(:,1); % extend time column
extratrop_landmo(last_ind_2000+1:last_ind_2006,1) = landusemo(last_ind_2000+1:last_ind_2006,1);
extratrop_landmo(last_ind_2000+1:last_ind_2006,2) = 0;

% other variables should just be loaded and passed