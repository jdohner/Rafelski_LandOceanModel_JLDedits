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

function [ff,landuse,landuseExtra] = getsourcesink_scale4(predict,start_year,end_year,year_vector);

if predict == 1

   
% load data thru 2009, 2006, 2000 - ff is monthly, lu is annual

% both vectors begin in 1800
landUse_2006 = csvread('landUse_1800-2006.csv');
landUseExtra_2000 = csvread('landUseExtra_1800-2000.csv');
load fossilFuel_1751-2009.mat;
ff_2009 = ff1; % already monthly resolution
% shortening ff vector to begin at start_year
ff1_start = find(ff_2009(:,1) == start_year);
%ff1_end = find(ff_2009(:,1) == end_year);
ff_2009 = ff_2009(ff1_start:end,:);

% load extended data thru 2016 - all annual
ff_2016 = csvread('fossilFuel_1959-2016.csv');
% this one isn't working for csvread, so reading in as text file
%landUse_2016 = csvread('landUse_1959-2016.csv'); 
fid = fopen('landUse_1959-2016.txt');
C = textscan(fid,'%f %f', 'delimiter','\t');
fclose(fid);
landUse_2016(:,1) = C{1};
landUse_2016(:,2) = C{2};


% interpolate to monthly

% landuse 2006 to monthly
month_2006 = 1800:(1/12):2006;
landUse_monthly2006 = interp1(landUse_2006(:,1),landUse_2006(:,12),month_2006);
landusemo_2006(:,1) = month_2006;
landusemo_2006(:,2) = landUse_monthly2006; % value in ppm

% extratrop landuse 2000 to monthly
month_2000 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
landUseExtra_monthly2000 = interp1(landUseExtra_2000(:,1),landUseExtra_2000(:,2),month_2000);
extratrop_landmo2000(:,1) = month_2000;
extratrop_landmo2000(:,2) = landUseExtra_monthly2000; % value in ppm

% fossil fuel 2009 already in monthly

% fossil fuel 2016 to monthly
month_2016 = 1959:(1/12):2016;
ff2016_monthly = (interp1(ff_2016(:,1),ff_2016(:,2),month_2016)).';
ffmo_2016(:,1) = month_2016;
ffmo_2016(:,2) = ff2016_monthly;

% land use 2016 to monthly
LU2016_monthly = (interp1(landUse_2016(:,1),landUse_2016(:,2),month_2016)).';
lumo_2016(:,1) = month_2016;
lumo_2016(:,2) = LU2016_monthly;



% extend (making full length vectors (1800-2016)) and patch

full_monthly = start_year:(1/12):end_year;
i_1959 = find(ff_2009(:,1) == 1959); % same for ff, landuse, extraLU
ff(:,1) = full_monthly;
ff(1:i_1959-1,2) = ff_2009(1:i_1959-1,2);
ff(i_1959:end,2) = ffmo_2016(:,2);

landuse(:,1) = full_monthly; 
landuse(1:i_1959-1,2) = landusemo_2006(1:i_1959-1,2);
landuse(i_1959:end,2) = lumo_2016(:,2);

% don't have updated extratropical land use data, so just extend from 2000
% onward using 0
landuseExtra(:,1) = full_monthly; 
landuseExtra(1:length(extratrop_landmo2000),2) = extratrop_landmo2000(:,2);
landuseExtra(length(extratrop_landmo2000)+1:end,2) = 0;



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

month_2006 = 1800:(1/12):2006;
landmonth_2006 = interp1(landnowppm_JLD(:,1),landnowppm_JLD(:,12),month_2006);

month_2000 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
extralandmonth_2000 = interp1(extratrop_landppm_JLD(:,1),extratrop_landppm_JLD(:,2),month_2000);

landusemo_2006(:,1) = month_2006;
landusemo_2006(:,2) = landmonth_2006; % value in ppm

extratrop_landmo2000(:,1) = month_2000;
extratrop_landmo2000(:,2) = extralandmonth_2000; % value in ppm

last_ind_2006 = length(month_2006);
last_ind_2000 = length(month_2000);

%Extend extratropical emissions by assuming emissions are zero
%extratrop_landmo(:,1) = landusemo(:,1); % extend time column
extratrop_landmo2000(last_ind_2000+1:last_ind_2006,1) = landusemo_2006(last_ind_2000+1:last_ind_2006,1);
extratrop_landmo2000(last_ind_2000+1:last_ind_2006,2) = 0;

% other variables should just be loaded and passed
% other variables should just be loaded and passed
end