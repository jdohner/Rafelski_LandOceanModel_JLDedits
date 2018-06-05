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

clear all

% load usa, china, russia data

% load western europe data

% total global emissions

% everything else

figure
plot(china_total(:,1),china_total(:,2),usa_total(:,1),usa_total(:,2), ...
    russia_total(:,1),russia_total(:,2));
title('total emissions by nation');
legend('china','usa','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
xlim([1800 2020])
grid

% Coal is a solid ff formed over millions of yrs by decay of land vegetation. 
figure
subplot(5,1,1)
plot(china_solid(:,1),china_solid(:,2),usa_solid(:,1),usa_solid(:,2), ...
    russia_solid(:,1),russia_solid(:,2));
title('solid (coal) emissions by nation');
legend('china','usa','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
xlim([1900 2020])
grid

% Oil is a liquid ff, formed from the remains of marine microorgs deposited on the sea floor.
subplot(5,1,2)
plot(china_liquid(:,1),china_liquid(:,2),usa_liquid(:,1),usa_liquid(:,2), ...
    russia_liquid(:,1),russia_liquid(:,2));
title('liquid (oil) emissions by nation');
legend('china','usa','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])

% Natural gas is a gaseous ff, relatively clean compared to coal and oil. 
subplot(5,1,3)
plot(china_gas(:,1),china_gas(:,2),usa_gas(:,1),usa_gas(:,2), ...
    russia_gas(:,1),russia_gas(:,2));
title('gas (natural gas) emissions by nation');
legend('china','usa','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])

% co2 emissions through calcination rxn to produce cement
subplot(5,1,4)
plot(china_cement(:,1),china_cement(:,2),usa_cement(:,1),usa_cement(:,2), ...
    russia_cement(:,1),russia_cement(:,2));
title('cement emissions by nation');
legend('china','usa','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])

% burning off flammable gas released by pressure relief valves in
% industrial plants
subplot(5,1,5)
plot(china_flaring(:,1),china_flaring(:,2),usa_flaring(:,1),usa_flaring(:,2), ...
    russia_flaring(:,1),russia_flaring(:,2));
title('flaring emissions by nation');
legend('china','usa','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
grid
xlim([1900 2020])
