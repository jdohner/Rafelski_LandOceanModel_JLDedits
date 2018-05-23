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
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match'));

% load Houghton LU data through 2016
% note: Houghton gave me 1850-2015 record, used 2016 Houghton value as
% listed in GCPv1.3
LUhough = csvread('HoughLU_perscomm_2016.csv');
luYear = LUhough(1,1):(1/ts):LUhough(end,1);
LUhoughmo_0 = (interp1(LUhough(:,1),LUhough(:,2),luYear)).';
LUhoughmo(:,1) = luYear;
LUhoughmo(:,2) = LUhoughmo_0*d1*d; %convert from TgC to PgC to ppm

% load Hansis LU data 1849-2016
% record used for GCP, not one from Hasis et. al 2015
LUhansis = csvread('Pongratz2016_GCP_meanPasture_peat.csv'); 
luYear2 = LUhansis(1,1):(1/ts):LUhansis(end,1);
LUhansismo_0 = (interp1(LUhansis(:,1),LUhansis(:,2),luYear2)).';
LUhansismo(:,1) = luYear2;
LUhansismo(:,2) = LUhansismo_0*d; %convert from PgC to ppm

% load GCB averaged historical LU data
% lu data, TgC/yr
LUhist = csvread('GCPv1.3_historicalLU2016.csv');
luYear3 = LUhist(1,1):(1/ts):LUhist(end,1);
LUhistmo_0 = (interp1(LUhist(:,1),LUhist(:,2),luYear3)).';
LUhistmo(:,1) = luYear3;
LUhistmo(:,2) = LUhistmo_0*d; %convert from PgC to ppm

% load LR high & low land use 
load LR_LU_LUex2016.mat; % land use records used in Rafelski 2009

% plot all LU cases over one another
figure
subplot(2,2,[1 2])
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2), ...
    LUhistmo(:,1),LUhistmo(:,2),LU(:,1),LU(:,2), LUex(:,1),LUex(:,2))
grid
title('Land Use Datasets', 'FontSize', 24)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 14)
xlabel('year', 'FontSize', 24)
ylabel('ppm', 'FontSize', 24)

% load residuals from driver using most recent ff data, all different lu
% cases 5 total - all for CHM-V with constant SST
load obsCO2_record.mat;
load updatedYear_vec.mat; % gives year2 vector
load updatedHough_outputvars.mat;
hough_co2 = atmcalc2;
hough_resid = obsCalcDiff;
load updatedhansis_outputvars.mat;
hansis_co2 = atmcalc2;
hansis_resid = obsCalcDiff;
load updatedGCPhist_outputvars.mat;
gcp_co2 = atmcalc2;
gcp_resid = obsCalcDiff;
load updatedLRLU_outputvars.mat;
LRLU_co2 = atmcalc2;
LRLU_resid = obsCalcDiff;
load updatedLRLUex_outputvars.mat;
LRLUex_co2 = atmcalc2;
LRLUex_resid = obsCalcDiff;

% plot modeled co2 records over one another
subplot(2,2,3)
plot(year2,CO2a(:,2),year2,hough_co2,year2,hansis_co2,year2,gcp_co2,year2,...
    LRLU_co2,year2,LRLUex_co2)
grid
title('Modeled and Observed CO_2 Records', 'FontSize', 24)
legend({'Observed CO2','Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 14)
xlabel('year', 'FontSize', 24)
ylabel('ppm', 'FontSize', 24)

% plot residuals over one another
subplot(2,2,4)
plot(year2,hough_resid(:,2),year2,hansis_resid(:,2),year2,gcp_resid(:,2),...
    year2,LRLU_resid(:,2),year2,LRLUex_resid(:,2))
grid
title('Observed - Modeled CO_2', 'FontSize', 24)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'FontSize', 14)
line([year2(1),year2(end)],[0,0],'linestyle','--');
xlabel('year', 'FontSize', 24)
ylabel('ppm', 'FontSize', 24)


% smooth residual of Houghton

% plot derivative