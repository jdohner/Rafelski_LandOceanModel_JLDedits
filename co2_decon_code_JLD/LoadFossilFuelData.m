% file LoadFossilFuelData.m
% author Lauren Rafelski, modified by Julia Dohner
% brief Loads fossil fuel emissions data from BP Statistical Review of 
% World Energy (BP, 2008)

% 3/23/09: Change input file to BP_extrap_CDIAC_data_2007.xls
% 1/10/11: Change input file to BP_extrap_CDIAC_data_2009.xls

% TODO: Does Marland et al. (2006) appear anywhere in this?
% Another option: give compact variable names, but be clear in defining
% them at the beginning, and comment a lot
% 9/16: last thing I did was line 38 (create new time vector)

function [ff1] = LoadFossilFuelData(timeStepPerYear)

%% load fossil fuel data

%Read in the data to be interpolated
fossilFuel=xlsread('BP_extrap_CDIAC_data_2009.xls');
fossilFuel_year = fossilFuel(:,1); %first column of .xls file is year
fossilFuel_emission = fossilFuel(:,2); %second column of .xls file is emission data

%Get rid of any NaN values
fossilFuel_emission(isnan(fossilFuel_year)) = [];
fossilFuel_year(isnan(fossilFuel_year)) = [];

%Compute cumulative flux
cumulativeFlux = 0;
fossilFuel_emission_cumulative = fossilFuel_emission;
yr1 = fossilFuel_year;
for i = 1:length(fossilFuel_year) %Loop through emissions matrix to calculate the cumulative FF emissions for
%each year
    cumulativeFlux = cumulativeFlux + fossilFuel_emission(i);
    fossilFuel_emission_cumulative(i) = cumulativeFlux;
    yr1(i) = fossilFuel_year(i)+1;  %add 1 because integral is valid at the end of the calendar year
end

%Create new time array
yr2 = fossilFuel_year(1):(1/timeStepPerYear):fossilFuel_year(end)+1;

%Do interpolation
fos2 = interp1(yr1,fossilFuel_emission_cumulative,yr2,'spline');

%Computed fluxes by taking discrete time derivative of integrated fluxes
fos = fos2;  %creates dummy arrays
yr = yr2;
fos(length(yr2)) = [];
yr(length(yr2)) = [];
for i = 1:(length(yr2)-1),
    fos(i) = timeStepPerYear.*(fos2(i+1)-fos2(i));
    yr(i) = yr2(i);
end

% convert to ppm

fosppm = fos/2.12;

ff1(:,2) = fosppm(1,:);
ff1(:,1) = yr(1,:);
