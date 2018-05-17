% getSourceSink5.m
%
% author: Lauren Rafelski, modified by Julia Dohner
% May 16, 2018
%
% gss version that uses most recent values for FF (Boden 2015, GCP) and 
% land use (Houghton 2015, personal comm.)

function [ff,LU] = getSourceSink5(year, ts);

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/necessary_data'));

d = 1/2.31; % gigaton to ppm conversion factor
d1 = 0.001; % teragram to petagram conversion factor
    
% ff data 1750-2015, GtC/yr
FF_2015 = csvread('FF_Boden_2015.csv');
ffYear = FF_2015(1,1):(1/ts):FF_2015(end,1);
FF_2015mo_0 = (interp1(FF_2015(:,1),FF_2015(:,2),ffYear)).';
FF_2015mo(:,1) = ffYear;
FF_2015mo(:,2) = FF_2015mo_0*d; %convert to pp

% lu data, TgC/yr
LU_2015 = csvread('Houghton_global_2015.csv');
luYear = LU_2015(1,1):(1/ts):LU_2015(end,1);
LU_2015mo_0 = (interp1(LU_2015(:,1),LU_2015(:,2),luYear)).';
LU_2015mo(:,1) = luYear;
LU_2015mo(:,2) = LU_2015mo_0*d1*d; %convert from TgC to PgC to ppm

% shorten datasets to match time frame of year vector
FF_start = find(FF_2015mo(:,1) == year(1));
FF_end = find(FF_2015mo(:,1) == year(end));
ff = FF_2015mo(FF_start:FF_end,:);
LU_start = find(LU_2015mo(:,1) == year(1));
LU_end = find(LU_2015mo(:,1) == year(end));
LU = LU_2015mo(LU_start:LU_end,:);

end