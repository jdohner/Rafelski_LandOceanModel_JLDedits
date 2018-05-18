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
    
% ff data 1750-2016, GtC/yr
FF_2016 = csvread('GCPv1.3_FF2016.csv');
ffYear = FF_2016(1,1):(1/ts):FF_2016(end,1);
FF_2016mo_0 = (interp1(FF_2016(:,1),FF_2016(:,2),ffYear)).';
FF_2016mo(:,1) = ffYear;
FF_2016mo(:,2) = FF_2016mo_0*d; %convert to pp

% lu data, TgC/yr
%LU_2016 = csvread('GCPv1.3_historicalLU2016.csv');
LU_2015 = csvread('Houghton_global_2015.csv');
luYear = LU_2016(1,1):(1/ts):LU_2016(end,1);
LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),luYear)).';
LU_2016mo(:,1) = luYear;
LU_2016mo(:,2) = LU_2016mo_0*d1*d; %convert from TgC to PgC to ppm

% shorten datasets to match time frame of year vector
FF_start = find(FF_2016mo(:,1) == year(1));
FF_end = find(FF_2016mo(:,1) == year(end));
ff = FF_2016mo(FF_start:FF_end,:);
LU_start = find(LU_2016mo(:,1) == year(1));
LU_end = find(LU_2016mo(:,1) == year(end));
LU = LU_2016mo(LU_start:LU_end,:);

end