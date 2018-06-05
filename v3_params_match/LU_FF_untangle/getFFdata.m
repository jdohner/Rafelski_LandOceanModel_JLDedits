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

WEstart = 1899; % constrained by china
WEend = 2013; % constrained by china

year_ann = WEstart:WEend;

china = csvread('china_emissions.csv');
china_total = [china(:,1) china(:,2)];
china_solid = [china(:,1) china(:,3)];
china_liquid = [china(:,1) china(:,4)];
china_gas = [china(:,1) china(:,5)];
china_cement = [china(:,1) china(:,6)];
china_flaring = [china(:,1) china(:,7)];

china_interp0 = (interp1(china_total(:,1),china_total(:,2),year_ann)).';
china_totalann(:,1) = year_ann;
china_totalann(:,2) = china_interp0;

usa = csvread('usa_emissions.csv');
usa_total = [usa(:,1) usa(:,2)];
usa_solid = [usa(:,1) usa(:,3)];
usa_liquid = [usa(:,1) usa(:,4)];
usa_gas = [usa(:,1) usa(:,5)];
usa_cement = [usa(:,1) usa(:,6)];
usa_flaring = [usa(:,1) usa(:,7)];

usa_interp0 = (interp1(usa_total(:,1),usa_total(:,2),year_ann)).';
usa_totalann(:,1) = year_ann;
usa_totalann(:,2) = usa_interp0;

russia = csvread('ussr_russia_emissions.csv');
russia_total = [russia(:,1) russia(:,2)];
russia_solid = [russia(:,1) russia(:,3)];
russia_liquid = [russia(:,1) russia(:,4)];
russia_gas = [russia(:,1) russia(:,5)];
russia_cement = [russia(:,1) russia(:,6)];
russia_flaring = [russia(:,1) russia(:,7)];

russia_interp0 = (interp1(russia_total(:,1),russia_total(:,2),year_ann)).';
russia_totalann(:,1) = year_ann;
russia_totalann(:,2) = russia_interp0;


%% Western Europe

belgium = csvread('belgium_emissions.csv');
belgium_total = [belgium(:,1) belgium(:,2)];
belgium_solid = [belgium(:,1) belgium(:,3)];
belgium_liquid = [belgium(:,1) belgium(:,4)];
belgium_gas = [belgium(:,1) belgium(:,5)];
belgium_cement = [belgium(:,1) belgium(:,6)];
belgium_flaring = [belgium(:,1) belgium(:,7)];

belgium_interp0 = (interp1(belgium_total(:,1),belgium_total(:,2),year_ann)).';
belgium_totalann(:,1) = year_ann;
belgium_totalann(:,2) = belgium_interp0;


france = csvread('france_emissions.csv');
france_total = [france(:,1) france(:,2)];
france_solid = [france(:,1) france(:,3)];
france_liquid = [france(:,1) france(:,4)];
france_gas = [france(:,1) france(:,5)];
france_cement = [france(:,1) france(:,6)];
france_flaring = [france(:,1) france(:,7)];

france_interp0 = (interp1(france_total(:,1),france_total(:,2),year_ann)).';
france_totalann(:,1) = year_ann;
france_totalann(:,2) = france_interp0;


% ireland = csvread('ireland_emissions.csv');
% ireland_total = [ireland(:,1) ireland(:,2)];
% ireland_solid = [ireland(:,1) ireland(:,3)];
% ireland_liquid = [ireland(:,1) ireland(:,4)];
% ireland_gas = [ireland(:,1) ireland(:,5)];
% ireland_cement = [ireland(:,1) ireland(:,6)];
% ireland_flaring = [ireland(:,1) ireland(:,7)];

% luxembourg = csvread('luxembourg_emissions.csv');
% luxembourg_total = [luxembourg(:,1) luxembourg(:,2)];
% luxembourg_solid = [luxembourg(:,1) luxembourg(:,3)];
% luxembourg_liquid = [luxembourg(:,1) luxembourg(:,4)];
% luxembourg_gas = [luxembourg(:,1) luxembourg(:,5)];
% luxembourg_cement = [luxembourg(:,1) luxembourg(:,6)];
% luxembourg_flaring = [luxembourg(:,1) luxembourg(:,7)];

netherlands = csvread('netherlands_emissions.csv');
netherlands_total = [netherlands(:,1) netherlands(:,2)];
netherlands_solid = [netherlands(:,1) netherlands(:,3)];
netherlands_liquid = [netherlands(:,1) netherlands(:,4)];
netherlands_gas = [netherlands(:,1) netherlands(:,5)];
netherlands_cement = [netherlands(:,1) netherlands(:,6)];
netherlands_flaring = [netherlands(:,1) netherlands(:,7)];

netherlands_interp0 = (interp1(netherlands_total(:,1),netherlands_total(:,2),year_ann)).';
netherlands_totalann(:,1) = year_ann;
netherlands_totalann(:,2) = netherlands_interp0;


uk = csvread('uk_emissions.csv');
uk_total = [uk(:,1) uk(:,2)];
uk_solid = [uk(:,1) uk(:,3)];
uk_liquid = [uk(:,1) uk(:,4)];
uk_gas = [uk(:,1) uk(:,5)];
uk_cement = [uk(:,1) uk(:,6)];
uk_flaring = [uk(:,1) uk(:,7)];

uk_interp0 = (interp1(uk_total(:,1),uk_total(:,2),year_ann)).';
uk_totalann(:,1) = year_ann;
uk_totalann(:,2) = uk_interp0;


italy = csvread('italy_emissions.csv');
italy_total = [italy(:,1) italy(:,2)];
italy_solid = [italy(:,1) italy(:,3)];
italy_liquid = [italy(:,1) italy(:,4)];
italy_gas = [italy(:,1) italy(:,5)];
italy_cement = [italy(:,1) italy(:,6)];
italy_flaring = [italy(:,1) italy(:,7)];

italy_interp0 = (interp1(italy_total(:,1),italy_total(:,2),year_ann)).';
italy_totalann(:,1) = year_ann;
italy_totalann(:,2) = italy_interp0;


germany = csvread('germany_full_emissions.csv');
germany_total = [germany(:,1) germany(:,2)];
germany_solid = [germany(:,1) germany(:,3)];
germany_liquid = [germany(:,1) germany(:,4)];
germany_gas = [germany(:,1) germany(:,5)];
germany_cement = [germany(:,1) germany(:,6)];
germany_flaring = [germany(:,1) germany(:,7)];

germany_interp0 = (interp1(germany_total(:,1),germany_total(:,2),year_ann)).';
germany_totalann(:,1) = year_ann;
germany_totalann(:,2) = germany_interp0;


spain = csvread('spain_emissions.csv');
spain_total = [spain(:,1) spain(:,2)];
spain_solid = [spain(:,1) spain(:,3)];
spain_liquid = [spain(:,1) spain(:,4)];
spain_gas = [spain(:,1) spain(:,5)];
spain_cement = [spain(:,1) spain(:,6)];
spain_flaring = [spain(:,1) spain(:,7)];

spain_interp0 = (interp1(spain_total(:,1),spain_total(:,2),year_ann)).';
spain_totalann(:,1) = year_ann;
spain_totalann(:,2) = spain_interp0;


% portugal = csvread('portugal_emissions.csv');
% portugal_total = [portugal(:,1) portugal(:,2)];
% portugal_solid = [portugal(:,1) portugal(:,3)];
% portugal_liquid = [portugal(:,1) portugal(:,4)];
% portugal_gas = [portugal(:,1) portugal(:,5)];
% portugal_cement = [portugal(:,1) portugal(:,6)];
% portugal_flaring = [portugal(:,1) portugal(:,7)];

norway = csvread('norway_emissions.csv');
norway_total = [norway(:,1) norway(:,2)];
norway_solid = [norway(:,1) norway(:,3)];
norway_liquid = [norway(:,1) norway(:,4)];
norway_gas = [norway(:,1) norway(:,5)];
norway_cement = [norway(:,1) norway(:,6)];
norway_flaring = [norway(:,1) norway(:,7)];

norway_interp0 = (interp1(norway_total(:,1),norway_total(:,2),year_ann)).';
norway_totalann(:,1) = year_ann;
norway_totalann(:,2) = norway_interp0;


sweden = csvread('sweden_emissions.csv');
sweden_total = [sweden(:,1) sweden(:,2)];
sweden_solid = [sweden(:,1) sweden(:,3)];
sweden_liquid = [sweden(:,1) sweden(:,4)];
sweden_gas = [sweden(:,1) sweden(:,5)];
sweden_cement = [sweden(:,1) sweden(:,6)];
sweden_flaring = [sweden(:,1) sweden(:,7)];

sweden_interp0 = (interp1(sweden_total(:,1),sweden_total(:,2),year_ann)).';
sweden_totalann(:,1) = year_ann;
sweden_totalann(:,2) = sweden_interp0;


%% total emissions - data processing

% million metric tons of C
% everything else in thousand metric tons
% 1 thousand metric ton C = 10e-6 PgC
d = 1000; % convert to thousand metric tons
global_data = csvread('global_total.csv');
global_total = [global_data(:,1) global_data(:,2).*d];
global_solid = [global_data(:,1) global_data(:,3).*d];
global_liquid = [global_data(:,1) global_data(:,4).*d];
global_gas = [global_data(:,1) global_data(:,5).*d];
global_cement = [global_data(:,1) global_data(:,6).*d];
global_flaring = [global_data(:,1) global_data(:,7).*d];



%% calculate Western Europe & Japan


belgium_cut = belgium_totalann((find(belgium_totalann(:,1) == WEstart)):(find(belgium_totalann(:,1) == WEend)),:);
france_cut = france_totalann((find(france_totalann(:,1) == WEstart)):(find(france_totalann(:,1) == WEend)),:);
%ireland_cut = ireland_total((find(ireland_total(:,1) == WEstart)):end,:);
italy_cut = italy_totalann((find(italy_totalann(:,1) == WEstart)):(find(italy_totalann(:,1) == WEend)),:);
netherlands_cut = netherlands_totalann((find(netherlands_totalann(:,1) == WEstart)):(find(netherlands_totalann(:,1) == WEend)),:);
norway_cut = norway_totalann((find(norway_totalann(:,1) == WEstart)):(find(norway_totalann(:,1) == WEend)),:);
sweden_cut = sweden_totalann((find(sweden_totalann(:,1) == WEstart)):(find(sweden_totalann(:,1) == WEend)),:);
uk_cut = uk_totalann((find(uk_totalann(:,1) == WEstart)):(find(uk_totalann(:,1) == WEend)),:);
germany_cut = germany_totalann((find(germany_totalann(:,1) == WEstart)):(find(germany_totalann(:,1) == WEend)),:);
spain_cut = spain_totalann((find(spain_totalann(:,1) == WEstart)):(find(spain_totalann(:,1) == WEend)),:);

japan = csvread('japan_emissions.csv');
japan_total = [japan(:,1) japan(:,2)];
japan_solid = [japan(:,1) japan(:,3)];
japan_liquid = [japan(:,1) japan(:,4)];
japan_gas = [japan(:,1) japan(:,5)];
japan_cement = [japan(:,1) japan(:,6)];
japan_flaring = [japan(:,1) japan(:,7)];

japan_interp0 = (interp1(japan_total(:,1),japan_total(:,2),year_ann)).';
japan_totalann(:,1) = year_ann;
japan_totalann(:,2) = japan_interp0;


japan_cut = japan_totalann((find(japan_totalann(:,1) == WEstart)):(find(japan_totalann(:,1) == WEend)),:);

WE_japan_total(:,1) = WEstart:WEend;
WE_japan_total(:,2) = belgium_cut(:,2)+france_cut(:,2)...
    + italy_cut(:,2) + netherlands_cut(:,2) + norway_cut(:,2)...
    + sweden_cut(:,2) + uk_cut(:,2) + japan_cut(:,2);

%% calculate "everything else"

usa_cut = usa_totalann((find(usa_totalann(:,1) == WEstart)):(find(usa_totalann(:,1) == WEend)),:);
china_cut = china_totalann((find(china_totalann(:,1) == WEstart)):(find(china_totalann(:,1) == WEend)),:);
russia_cut = russia_totalann((find(russia_totalann(:,1) == WEstart)):(find(russia_totalann(:,1) == WEend)),:);

everythingElse(:,1) = WEstart:WEend;
global_cut = global_total((find(global_total(:,1) == WEstart)):(find(global_total(:,1) == WEend)),:);
everythingElse(:,2) = global_cut(:,2) - WE_japan_total(:,2) - usa_cut(:,2)...
    -china_cut(:,2)-russia_cut(:,2);

%% replace missing values with nan's

save('ffdata_total','china_totalann','usa_totalann','russia_totalann',...
    'global_total','WE_japan_total','everythingElse')

