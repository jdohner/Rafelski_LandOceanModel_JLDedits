% LU_interp.mat
%
% script to interploate LU records to monthly
%
% author: julia dohner

clear all

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;
year = (start_year:(1/ts):end_year)';
d = 2.31; % ppm to PgC conversion factor (formerly 1/2.31 opp direction)
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
LUhoughmo(:,2) = LUhoughmo_0*d1; %convert from TgC to PgC

% load Hansis LU data 1849-2016
% record used for GCP, not one from Hasis et. al 2015
LUhansis = csvread('Pongratz2016_GCP_meanPasture_peat.csv'); 
luYear2 = LUhansis(1,1):(1/ts):LUhansis(end,1);
LUhansismo_0 = (interp1(LUhansis(:,1),LUhansis(:,2),luYear2)).';
LUhansismo(:,1) = luYear2;
LUhansismo(:,2) = LUhansismo_0; 

% load GCB averaged historical LU data
% lu data, TgC/yr
LUhist = csvread('GCPv1.3_historicalLU2016.csv');
luYear3 = LUhist(1,1):(1/ts):LUhist(end,1);
LUgcpmo_0 = (interp1(LUhist(:,1),LUhist(:,2),luYear3)).';
LUgcpmo(:,1) = luYear3;
LUgcpmo(:,2) = LUgcpmo_0; 

% load LR high & low land use 
load LR_LU_LUex2016.mat; % land use records used in Rafelski 2009
LUhough03mo(:,2) = LU(:,2)*d; % convert from ppm to PgC
%LUex(:,2) = LUex(:,2)*d; % convert from ppm to PgC

save('LU_records_monthly','LUhansismo','LUgcpmo','LUhoughmo','LUhough03mo','year');