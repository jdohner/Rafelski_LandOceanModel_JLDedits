% plot_national_ff.mat
%
% script to plot different types of emissions from different countries
%
% julia dohner
% may 15, 2018
%
% 						
% Data Source: Tom Boden and Bob Andres (Oak Ridge National Laboratory); 
% Gregg Marland (Appalachian State University)									
% DOI: 10.3334/CDIAC/00001_V2017	


% format of columns being read in from .csv files:
%Nation	
%Year	
%Total CO2 emissions from fossil-fuels and cement production (thousand metric tons of C)	
%Emissions from solid fuel consumption	
%Emissions from liquid fuel consumption	
%Emissions from gas fuel consumption	
%Emissions from cement production	
%Emissions from gas flaring	Per capita CO2 emissions (metric tons of carbon)	
%Emissions from bunker fuels (not included in the totals)

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle/nationalFFdata'));
% top 10: 1. china, 2. usa, 3. india, 4. russia, 5. japan, 6. germany

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

india = csvread('india_emissions.csv');
india_total = [india(:,1) india(:,2)];
india_solid = [india(:,1) india(:,3)];
india_liquid = [india(:,1) india(:,4)];
india_gas = [india(:,1) india(:,5)];
india_cement = [india(:,1) india(:,6)];
india_flaring = [india(:,1) india(:,7)];

russia = csvread('ussr_russia_emissions.csv');
russia_total = [russia(:,1) russia(:,2)];
russia_solid = [russia(:,1) russia(:,2)];
russia_liquid = [russia(:,1) russia(:,2)];
russia_gas = [russia(:,1) russia(:,2)];
russia_cement = [russia(:,1) russia(:,2)];
russia_flaring = [russia(:,1) russia(:,2)];

figure
plot(china_total(:,1),china_total(:,2),usa_total(:,1),usa_total(:,2), ...
    india_total(:,1),india_total(:,2),russia_total(:,1),russia_total(:,2));
title('total emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
xlim([1800 2020])
grid

% Coal is a solid ff formed over millions of yrs by decay of land vegetation. 
figure
subplot(5,1,1)
plot(china_solid(:,1),china_solid(:,2),usa_solid(:,1),usa_solid(:,2), ...
    india_solid(:,1),india_solid(:,2),russia_solid(:,1),russia_solid(:,2));
title('solid (coal) emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
xlim([1900 2020])
grid

% Oil is a liquid ff, formed from the remains of marine microorgs deposited on the sea floor.
subplot(5,1,2)
plot(china_liquid(:,1),china_liquid(:,2),usa_liquid(:,1),usa_liquid(:,2), ...
    india_liquid(:,1),india_liquid(:,2),russia_liquid(:,1),russia_liquid(:,2));
title('liquid (oil) emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])

% Natural gas is a gaseous ff, relatively clean compared to coal and oil. 
subplot(5,1,3)
plot(china_gas(:,1),china_gas(:,2),usa_gas(:,1),usa_gas(:,2), ...
    india_gas(:,1),india_gas(:,2),russia_gas(:,1),russia_gas(:,2));
title('gas (natural gas) emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])

% co2 emissions through calcination rxn to produce cement
subplot(5,1,4)
plot(china_cement(:,1),china_cement(:,2),usa_cement(:,1),usa_cement(:,2), ...
    india_cement(:,1),india_cement(:,2),russia_cement(:,1),russia_cement(:,2));
title('cement emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])

% burning off flammable gas released by pressure relief valves in
% industrial plants
subplot(5,1,5)
plot(china_flaring(:,1),china_flaring(:,2),usa_flaring(:,1),usa_flaring(:,2), ...
    india_flaring(:,1),india_flaring(:,2),russia_flaring(:,1),russia_flaring(:,2));
title('flaring emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])
