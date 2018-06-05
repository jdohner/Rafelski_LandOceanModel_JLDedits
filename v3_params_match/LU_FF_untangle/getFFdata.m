% getFFdata.m
%
% script for prepping western europe ff data
%
% julia dohner
% june 4, 2018


clear all

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle/nationalFFdata'));

% units in thousand metric tons
% 1 thousand metric ton C = 10e-6 PgC


china = csvread('china_emissions.csv');
china_total = [china(:,1) china(:,2)];
china_solid = [china(:,1) china(:,3)];
china_liquid = [china(:,1) china(:,4)];
china_gas = [china(:,1) china(:,5)];
china_cement = [china(:,1) china(:,6)];
china_flaring = [china(:,1) china(:,7)];

usa = csvread('usa_emissions.csv');
usa_total = [usa(:,1) usa(:,2)];
usa_solid = [usa(:,1) usa(:,3)];
usa_liquid = [usa(:,1) usa(:,4)];
usa_gas = [usa(:,1) usa(:,5)];
usa_cement = [usa(:,1) usa(:,6)];
usa_flaring = [usa(:,1) usa(:,7)];

russia = csvread('ussr_russia_emissions.csv');
russia_total = [russia(:,1) russia(:,2)];
russia_solid = [russia(:,1) russia(:,3)];
russia_liquid = [russia(:,1) russia(:,4)];
russia_gas = [russia(:,1) russia(:,5)];
russia_cement = [russia(:,1) russia(:,6)];
russia_flaring = [russia(:,1) russia(:,7)];

%% Western Europe

belgium = csvread('belgium_emissions.csv');
belgium_total = [belgium(:,1) belgium(:,2)];
belgium_solid = [belgium(:,1) belgium(:,3)];
belgium_liquid = [belgium(:,1) belgium(:,4)];
belgium_gas = [belgium(:,1) belgium(:,5)];
belgium_cement = [belgium(:,1) belgium(:,6)];
belgium_flaring = [belgium(:,1) belgium(:,7)];


france = csvread('france_emissions.csv');
france_total = [france(:,1) france(:,2)];
france_solid = [france(:,1) france(:,3)];
france_liquid = [france(:,1) france(:,4)];
france_gas = [france(:,1) france(:,5)];
france_cement = [france(:,1) france(:,6)];
france_flaring = [france(:,1) france(:,7)];

ireland = csvread('ireland_emissions.csv');
ireland_total = [ireland(:,1) ireland(:,2)];
ireland_solid = [ireland(:,1) ireland(:,3)];
ireland_liquid = [ireland(:,1) ireland(:,4)];
ireland_gas = [ireland(:,1) ireland(:,5)];
ireland_cement = [ireland(:,1) ireland(:,6)];
ireland_flaring = [ireland(:,1) ireland(:,7)];

luxembourg = csvread('luxembourg_emissions.csv');
luxembourg_total = [luxembourg(:,1) luxembourg(:,2)];
luxembourg_solid = [luxembourg(:,1) luxembourg(:,3)];
luxembourg_liquid = [luxembourg(:,1) luxembourg(:,4)];
luxembourg_gas = [luxembourg(:,1) luxembourg(:,5)];
luxembourg_cement = [luxembourg(:,1) luxembourg(:,6)];
luxembourg_flaring = [luxembourg(:,1) luxembourg(:,7)];

netherlands = csvread('netherlands_emissions.csv');
netherlands_total = [netherlands(:,1) netherlands(:,2)];
netherlands_solid = [netherlands(:,1) netherlands(:,3)];
netherlands_liquid = [netherlands(:,1) netherlands(:,4)];
netherlands_gas = [netherlands(:,1) netherlands(:,5)];
netherlands_cement = [netherlands(:,1) netherlands(:,6)];
netherlands_flaring = [netherlands(:,1) netherlands(:,7)];

uk = csvread('uk_emissions.csv');
uk_total = [uk(:,1) uk(:,2)];
uk_solid = [uk(:,1) uk(:,3)];
uk_liquid = [uk(:,1) uk(:,4)];
uk_gas = [uk(:,1) uk(:,5)];
uk_cement = [uk(:,1) uk(:,6)];
uk_flaring = [uk(:,1) uk(:,7)];

italy = csvread('italy_emissions.csv');
italy_total = [italy(:,1) italy(:,2)];
italy_solid = [italy(:,1) italy(:,3)];
italy_liquid = [italy(:,1) italy(:,4)];
italy_gas = [italy(:,1) italy(:,5)];
italy_cement = [italy(:,1) italy(:,6)];
italy_flaring = [italy(:,1) italy(:,7)];

germany = csvread('germany_full_emissions.csv');
germany_total = [germany(:,1) germany(:,2)];
germany_solid = [germany(:,1) germany(:,3)];
germany_liquid = [germany(:,1) germany(:,4)];
germany_gas = [germany(:,1) germany(:,5)];
germany_cement = [germany(:,1) germany(:,6)];
germany_flaring = [germany(:,1) germany(:,7)];

spain = csvread('spain_emissions.csv');
spain_total = [spain(:,1) spain(:,2)];
spain_solid = [spain(:,1) spain(:,3)];
spain_liquid = [spain(:,1) spain(:,4)];
spain_gas = [spain(:,1) spain(:,5)];
spain_cement = [spain(:,1) spain(:,6)];
spain_flaring = [spain(:,1) spain(:,7)];

portugal = csvread('portugal_emissions.csv');
portugal_total = [portugal(:,1) portugal(:,2)];
portugal_solid = [portugal(:,1) portugal(:,3)];
portugal_liquid = [portugal(:,1) portugal(:,4)];
portugal_gas = [portugal(:,1) portugal(:,5)];
portugal_cement = [portugal(:,1) portugal(:,6)];
portugal_flaring = [portugal(:,1) portugal(:,7)];

norway = csvread('norway_emissions.csv');
norway_total = [norway(:,1) norway(:,2)];
norway_solid = [norway(:,1) norway(:,3)];
norway_liquid = [norway(:,1) norway(:,4)];
norway_gas = [norway(:,1) norway(:,5)];
norway_cement = [norway(:,1) norway(:,6)];
norway_flaring = [norway(:,1) norway(:,7)];

sweden = csvread('sweden_emissions.csv');
sweden_total = [sweden(:,1) sweden(:,2)];
sweden_solid = [sweden(:,1) sweden(:,3)];
sweden_liquid = [sweden(:,1) sweden(:,4)];
sweden_gas = [sweden(:,1) sweden(:,5)];
sweden_cement = [sweden(:,1) sweden(:,6)];
sweden_flaring = [sweden(:,1) sweden(:,7)];

%% total emissions - data processing

global_data = csvread('global_total.csv');
global_total = [global_data(:,1) global_data(:,2)];
global_solid = [global_data(:,1) global_data(:,3)];
global_liquid = [global_data(:,1) global_data(:,4)];
global_gas = [global_data(:,1) global_data(:,5)];
global_cement = [global_data(:,1) global_data(:,6)];
global_flaring = [global_data(:,1) global_data(:,7)];

%% calculate Western Europe & Japan



%% calculate "everything else"



%% replace missing values with nan's

save('ffdata_total','china_total','usa_total','russia_total',...
    'global_total','WE_japan_total','everythingElse')

