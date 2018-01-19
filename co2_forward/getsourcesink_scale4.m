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

function [ff,landusemo,extratrop_landmo] = getsourcesink_scale4(predict,start_year,end_year,year_vector);

if predict == 1

% LR data to 2000, 2006
load fossilFuel_1751-2009.mat;
ff_2009 = ff1;
landUse_2006 = csvread('landUse_1800-2006.csv');
landUseExtra_2000 = csvread('landUseExtra_1800-2000.csv');

% extended data
ff_2016 = csvread('fossilFuel_1959-2016.csv');
% interpolate to monthly
month_2016 = 1959:(1/12):2016;
ff2016_monthly = (interp1(ff_2016(:,1),ff_2016(:,2),month_2016)).';
ffmo_2016(:,1) = month_2016;
ffmo_2016(:,2) = ff2016_monthly;

% this one isn't working for csvread, so reading in as text file
%landUse_2016 = csvread('landUse_1959-2016.csv'); 
fid = fopen('landUse_1959-2016.txt');
C = textscan(fid,'%f %f', 'delimiter','\t');
fclose(fid);
landUse_2016(:,1) = C{1};
landUse_2016(:,2) = C{2};
LU2016_monthly = (interp1(landUse_2016(:,1),landUse_2016(:,2),month_2016)).';
lumo_2016(:,1) = month_2016;
lumo_2016(:,2) = LU2016_monthly;

% shortening ff vector to begin at start_year
ff1_start = find(ff_2009(:,1) == start_year);
ff1_end = find(ff_2009(:,1) == end_year);
ff = ff_2009(ff1_start:ff1_end,:);

% interpolate land use changes to monthly resolution
month = 1800:(1/12):2006;
landUse_monthly = interp1(landUse_2006(:,1),landUse_2006(:,12),month);

month2 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
landUseExtra_monthly = interp1(landUseExtra_2000(:,1),landUseExtra_2000(:,2),month2);

landusemo(:,1) = month;
landusemo(:,2) = landUse_monthly; % value in ppm

extratrop_landmo(:,1) = month2;
extratrop_landmo(:,2) = landUseExtra_monthly; % value in ppm

last_ind_2006 = length(month);
last_ind_2000 = length(month2);

%% note: LR ff and LU only go to 2009, 2006, so need to make vectors longer
% to avoid dimension mismatch

% replace 1959 onwards for fossil fuel, land use, extra trop land use
%ff = zeros(length(year_vector),2);
i = find(ff_2009(:,1) == 1959);
%ff(1:i-1) = ff_2009(1:i-1);
ff(i:length(ff)) = ff_2016(:);


else % if predict == 0 (diagnostic case, not forward projecting)
    

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

load joos_hilda_2011.mat
landnowppm_JLD = csvread('landnowppm_JLD.csv');
extratrop_landppm_JLD = csvread('extratrop_landppm_JLD.csv');


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
% other variables should just be loaded and passed
end