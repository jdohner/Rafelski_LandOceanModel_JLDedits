% file OceanUptakeModel_Driver.m
% formerly "jooshildascale_annotate2.m"
% 
% author Lauren Rafelski, modified by Julia Dohner
% 
% note Be sure to run defaults.m in outer folder before running this code
% 
% brief OceanUptakeModel_Driver.m is the main code that drives the subroutines.
% Together, these codes will calculate the ocean uptake of carbon dioxide
% using a pulse-response function from Joos et al. (1996). Then, the codes
% will calculate the residual land flux based on CO2 emissions, change in
% atmospheric CO2 concentration, and the ocean uptake.
% 
% Changes by LR:
% January 10, 2011: extended data to 2010
%
%
%TODO: compare runtime of this code to original code
%TODO: rename file as OceanUptake.m

clear all;

timeStepPerYear = 12; % number of data points/year
start_year = 1800;
end_year = 2010;
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
% Note: c converts air-sea flux from units of [ppm m^-2 yr^-1] into units of [umol
% m^-2 yr^-1], and converts m^3 into kg seawater (pg 402)


% Note: timeStepPerYear used to be called ts (NOT dt)
% fewer outputs taken in following line than provided by function
% get atmospheric CO2 record
[dtdelpCO2a,dpCO2a,year,dt] = GetIncrementMergedCO2(timeStepPerYear,start_year,end_year); 

% get fossil fuel emissions
[fossilFuelData] = LoadFossilFuelData(timeStepPerYear); 

% allocating space for t and r matrices before for loop
t = NaN(length(year), 1); % tracer concentration (Joos 1996 A.1. pg 415)
r = NaN(length(year), 2); % pulse response function (Joos 1996 A.1. pg 415)

% Response function to calculate ocean uptake
for n = 1:length(year)
     t(n,1) = year(n) - year(1);
     r(n,1) = t(n,1);
    %Calculate response function based on HILDA equation in Joos 1996
     if t(n,1) == 0
         r(n,2) = 1;
     elseif t(n,1) <= 2
         % A.2.2. HILDA model for 0<t<=2 yr (page 415)
         % (1/0.95873) multiplied in front
         % TODO: why multiplied by this?
         r(n,2)= (1/0.95873)*(0.12935+0.21898*exp(-t(n,1)/0.034569)+0.17003*exp(-t(n,1)/0.26936)...
             +0.24071*exp(-t(n,1)/0.96083)+0.24093*exp(-t(n,1)/4.9792));
     else
         % Joos (1996) A.2.2. HILDA model for 2 yr<t equation (page 415)
         % (1/0.95873) multiplied in front
         r(n,2) = (1/0.95873)*(0.022936+0.24278*exp(-t(n,1)/1.2679)+0.13963*exp(-t(n,1)/5.2528)...
                +0.089318*exp(-t(n,1)/18.601)+0.037820*exp(-t(n,1)/68.736)...
                +0.035549*exp(-t(n,1)/232.3));
     end
end

%% Calculate ocean uptake

[airSeaFlux,dpCO2s] = OceanPulseResponse(year,dpCO2a,c,h,kg,T,Aoc,r,dt); 

%% Calculate land flux using fossil fuel sources and ocean sink (ocean flux*area of ocean)

%allocate space for landFlux matrix
landFlux = NaN(length(year)-7, 2);

%TODO: finish commenting this
for p = 1:(length(year)-7) %1:(length(year)-1) %TODO: was this commented out by Lauren?
    q = find(fossilFuelData(:,1) == year(1,p));
    landFlux(p,1) = year(p);
    landFlux(p,2) = dtdelpCO2a(p+((1800-1640)*timeStepPerYear),2) - fossilFuelData(q,2) + airSeaFlux(p,2)*Aoc; 
end
% calculating LU + B (total land flux) in loop above -RFK

x = 0*year;

figure
plot(fossilFuelData(:,1),fossilFuelData(:,2),'-k',dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r',airSeaFlux(:,1),-Aoc*airSeaFlux(:,2),'-b',landFlux(:,1),landFlux(:,2),'-g',year(1,:),x,'--k')
axis([1800 2010 -10 10])
legend('fossil fuel','atmosphere','ocean','land','Location','SouthWest')
title('Sources and sinks from Joos response function ')
xlabel('Year ')
ylabel('ppm/year  Positive = source, negative = sink ')