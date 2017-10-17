% file GlobalCO2_Forward_Driver.m
% 
% author Julia Dohner, borrowing from Lauren Rafelski
% 
% note Be sure to run defaults.m in outer folder before running this code
% 
% brief 
% 
%
% subroutines: 
% HILDA response model r (calculates HILDA ocean uptake - A.2.2 in Joos 96)
% MLPulseResponse (calculates mixed layer uptake - Eq. 3, 6b in Joos 1996)
% biobox (Calculates co2 uptake of fast land box model (driven by co2
% fertilization))

clear all

% give access to data files in co2_forward_data folder
addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));

LU = 1; %1 = high land use scenario; 2 = low land use scenario

nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?

filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered


load land_temp.mat % land temperature records

load npp_T.mat % NPP-weighted temperature record

load landwt_T_2011.mat % land temperature anomaly

% load CO2 sources and sinks
% landusemo: land use emissions, in ppm/year, interpolated to monthly
% resolution
% ff1: fossil fuel emissions
% fas: ocean flux per m^2
% Aoc: ocean surface area
% extratrop_landmo: extratropical land use emissions

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3; 

 clear year start_year end_year ts

ts = 12; % timesteps per year
start_year = 1850;
end_year = 2009+(7/12); 

beta = [0.5;2]; % initial guesses for model fit

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year); 

% Extend land use record by making recent emissions equal to last
% record
landusemo(1874:1916,1) = year(1874:1916);
landusemo(1874:1916,2) = landusemo(1873,2);
% % 
%Extend extratropical emissions by assuming emissions are zero
extratrop_landmo(1802:1916,1) = landusemo(1802:1916,1);
extratrop_landmo(1802:1916,2) = 0;

%% Calculate residual land uptake
% run to 8/2009
% using high land use emissions
residual(:,1) = year(1,1:1916);
residual(:,2) = dtdelpCO2a(2521:4436,2) - ff1(1189:3104,2)....
+ Aoc*fas(601:2516,2) - landusemo(1:1916,2);

% using extratropical emissions only
residual2(:,1) = year(1,1:1916);
residual2(:,2) = dtdelpCO2a(2521:4436,2) - ff1(1189:3104,2)....
+ Aoc*fas(601:2516,2) - extratrop_landmo(1:1916,2);

