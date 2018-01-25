% Get relevant land use and Joos model results
%
% Updates
%
% 4/13/07 - Lauren Elmegreen - use land use change 
% emissions that are extrapolated to the present by keeping values 
% constant at the last value, instead of increasing by 1.4%
% 10/25/12 - Lauren Rafelski - added annotations
%
% This file is for the record between 1800 and 2006+10/12

function [ff1,landusemo,extratrop_landmo] = getsourcesink_scale5(startYr,endYr,dt);

load joos_hilda_2011.mat % contains year, start_year, end_year
landnowppm_JLD = csvread('landnowppm_JLD.csv');
extratrop_landppm_JLD = csvread('extratrop_landppm_JLD.csv');

year = startYr:dt:endYr;

% shortening ff vector to begin at start_year
ff1_start = find(ff1(:,1) == 1800);
ff1_end = find(ff1(:,1) == 2009+(7/12));
ff1 = ff1(ff1_start:ff1_end,:);

% interpolate land use changes to monthly resolution

month_2006 = 1800:(1/12):2006;
landmonth = interp1(landnowppm_JLD(:,1),landnowppm_JLD(:,12),month_2006);
landusemo(:,1) = month_2006;
landusemo(:,2) = landmonth; % value in ppm


month_2000 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
extralandmonth = interp1(extratrop_landppm_JLD(:,1),extratrop_landppm_JLD(:,2),month_2000);
extratrop_landmo(:,1) = month_2000;
extratrop_landmo(:,2) = extralandmonth; % value in ppm

% make new vectors through end of year

i_2000 = length(month_2000);
i_endYr = length(year);
% 
% % extend to end of year vector
% extratrop_landmo(i_2000+1:i_endYr,1) = year(i_2000+1:i_endYr,1);
% extratrop_landmo(i_2000+1:i_endYr,2) = 0;
% 
% %Extend extratropical emissions by assuming emissions are zero
% %extratrop_landmo(:,1) = landusemo(:,1); % extend time column
% extratrop_landmo(i_2000+1:i_endYr,1) = landusemo(i_2000+1:last_ind_2006,1);
% extratrop_landmo(last_ind_2000+1:last_ind_2006,2) = 0;
% 
% % other variables should just be loaded and passed


%% 

i_2006 = length(month_2006);
LU(:,1) = year;
LU(1:i_2006,2) = landusemo(1:i_2006,2);
LU(i_2006:end,2) = landusemo(i_2006,2);

i_2000 = length(month_2000);
LUex(:,1) = year;
LUex(1:i_2000,2) = extratrop_landmo(1:i_2000,2);


% 
% LU(:,1) = year; 
% LU(1:i_1959-1,2) = LU_2006mo(1:i_1959-1,2);
%     LU(i_1959:end,2) = LU_2016mo(:,2);
    
    
    
    
    % extend (making full length vectors (1800-2016)) and patch

    % Houghton record begins in 1959, combine two records at 1959
    i_1959 = find(FF_2009(:,1) == 1959); % same for ff, landuse, extraLU
    ff(:,1) = year;
    ff(1:i_1959-1,2) = FF_2009(1:i_1959-1,2);
    ff(i_1959:end,2) = FF_2016mo(:,2);

    LU(:,1) = year; 
    LU(1:i_1959-1,2) = LU_2006mo(1:i_1959-1,2);
    LU(i_1959:end,2) = LU_2016mo(:,2);

    % don't have updated extratropical land use data, so just extend from 2000
    % onward using 0
    LUex(:,1) = year; 
    LUex(1:length(LUex_2000mo),2) = LUex_2000mo(:,2);
    LUex(length(LUex_2000mo)+1:end,2) = 0;




%%%%
% getsourcesink outputs end at 2006
% Extend land use record by making recent emissions equal to last
% record
j = length(landuse);
k = length(year);
extendYears = year(1,j+1:k);
%landusemo(:,1) = [landusemo(:,1), extendYears];

landuse(j+1:k,1) = year(1,j+1:k);
landuse(j+1:k,2) = landuse(j,2);

%Extend extratropical emissions by assuming emissions are zero
landuseExtra(j+1:k,1) = landuse(j+1:k,1);
landuseExtra(j+1:k,2) = 0;