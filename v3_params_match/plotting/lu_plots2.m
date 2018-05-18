% lu_plots2.m
%
% script for plotting different land use data
%
% author: julia dohner
% may 4, 2018
%
% updated script using only full records (no appended)

clear all

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;
year2 = (start_year:(1/ts):end_year)';
d = 1/2.31; % gigaton to ppm conversion factor
d1 = 0.001; % teragram to petagram conversion factor


addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/necessary_data'));

% load Houghton LU data through 2016
% note: Houghton gave me 1850-2015 record, used 2016 Houghton value as
% listed in GCPv1.3
LUhough = csvread('HoughLU_perscomm_2016.csv');
luYear = LUhough(1,1):(1/ts):LUhough(end,1);
LUhoughmo_0 = (interp1(LUhough(:,1),LUhough(:,2),luYear)).';
LUhoughmo(:,1) = luYear;
LUhoughmo(:,2) = LUhoughmo_0*d1*d; %convert from TgC to PgC to ppm

% load GCB averaged historical LU data
% lu data, TgC/yr
LUhist = csvread('GCPv1.3_historicalLU2016.csv');
luYear = LUhist(1,1):(1/ts):LUhist(end,1);
LUhistmo_0 = (interp1(LUhist(:,1),LUhist(:,2),luYear)).';
LUhistmo(:,1) = luYear;
LUhistmo(:,2) = LUhistmo_0*d; %convert from TgC to PgC to ppm

% load LR high & low land use 
load LR_LU_LUex2016.mat; % land use records used in Rafelski 2009

% plot all LU cases over one another
figure
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhistmo(:,1),LUhistmo(:,2),LU(:,1),LU(:,2),...
    LUex(:,1),LUex(:,2))
legend('houghton','GCP','LR high','LR low')


% load residuals from driver using most recent ff data, all different lu
% cases

% plot residuals over one another

% smooth residual of Houghton

% plot derivative