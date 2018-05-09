% getSourceSink4.m
%
% author: Lauren Rafelski, modified by Julia Dohner
% January 23, 2018
%
% gss version that tries out different LU cases 
%
% LU in ppm/year
% ff in ppm/year
%

% which LU record to use before 1959????

function [ff,LU] = getSourceSink4(year2, ts);

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_data'));
d = 1/2.31; % gigaton to ppm conversion factor

year3 = (year2(1,1):(1/ts):2016)'; % making full data vector thru to 2016
    
% load data thru 2009, 2006, 2000 - ff is monthly, lu is annual



% both vectors begin in 1800
% LU_2006 = csvread('landUse_1800-2006.csv');
% LUex_2000 = csvread('landUseExtra_1800-2000.csv');
load fossilFuel_1751-2009.mat;
FF_2009 = ff1; % already monthly resolution
% shortening ff vector to begin at start_year (have data back thru 1700)
FF_start = find(FF_2009(:,1) == year2(1));
FF_2009 = FF_2009(FF_start:end,:);


% historical LU 1800-1960
LU_early = csvread('Historical_LU.csv');

% interpolate to monthly
earlyStart = LU_early(1,1);
earlyEnd = LU_early(end,1);
LU_earlymo(:,1) = (earlyStart:(1/12):earlyEnd)';
LU_earlymo(:,2) = (interp1(LU_early(:,1),LU_early(:,2),LU_earlymo(:,1)))*d;
LU_start = find(LU_earlymo(:,1) == year2(1)); % start at 1850
LU_early = LU_earlymo(LU_start:end,:); % end at 1960

%% extend vectors to 2016

% load extended data thru 2016 - all annual
FF_2016 = csvread('fossilFuel_1959-2016.csv'); % in gigatons/year
% 1 ppm CO2 = 2.31 gton CO2


% this one isn't working for csvread, so reading in as text file
%landUse_2016 = csvread('landUse_1959-2016.csv'); 
% fid = fopen('landUse_1959-2016.txt');
% C = textscan(fid,'%f %f', 'delimiter','\t');
% fclose(fid);
% LU_2016(:,1) = C{1};
% LU_2016(:,2) = C{2};

%GCB_LU.csv
%BLUE_LU.csv
%CABLE_LU.csv
%CLASS-CTEM_LU.csv
%Houghton_LU.csv
%LR_LU.csv
%LR_LUex.csv
LU_2016 = csvread('GCB_LU.csv');

% fossil fuel 2016 to monthly (arrives as 1959-2016)
month_2016 = 1959:(1/12):2016;
FF_2016mo_0 = (interp1(FF_2016(:,1),FF_2016(:,2),month_2016)).';
FF_2016mo(:,1) = month_2016;
FF_2016mo(:,2) = FF_2016mo_0*d; %convert to ppm

% land use 2016 to monthly
LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),month_2016)).';
LU_2016mo(:,1) = month_2016;
LU_2016mo(:,2) = LU_2016mo_0*d;

% extend (making full length vectors (1800-2016)) and patch

clear ff; % get rid of dimension mismatch since looks like below pulls
% from original data source (FF_2009) anyways
% Houghton record begins in 1959, combine two records at 1959
%i_2009 = find(FF_2009(:,1) == 2009); % same for ff, landuse, extraLU
i_2009 = length(FF_2009);
enddate_2009 = FF_2009(end,1);
startdate_2016 = find(FF_2016mo(:,1) == enddate_2009);
ff(:,1) = year3;
%j_2009 = find(FF_2009(:,1) == 2009);
ff(1:i_2009,2) = FF_2009(:,2);
ff(i_2009+1:end,2) = FF_2016mo(startdate_2016+1:end,2);

i_LU = find(LU_early(:,1) == LU_2016mo(1,1));
LU_extend = [LU_early(1:i_LU-1,:) ; LU_2016mo(:,:)];



%% shorten datasets to match time frame of year vector

FF_start = find(ff(:,1) == year2(1));
FF_end = find(ff(:,1) == year2(end));
ff = ff(FF_start:FF_end,:);
LU_start = find(LU_extend(:,1) == year2(1));
LU_end = find(LU_extend(:,1) == year2(end));
LU = LU_extend(LU_start:LU_end,:);
% LUex_end = find(LUex(:,1) == year2(end));
% LUex = LUex(1:LUex_end,:);