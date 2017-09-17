% file LoadFossilFuelData.m
% author Lauren Rafelski, modified by Julia Dohner
% brief Loads fossil fuel emissions data from BP Statistical Review of 
% World Energy (BP, 2008)

% data updates
% 3/23/09: Change input file to BP_extrap_CDIAC_data_2007.xls
% 1/10/11: Change input file to BP_extrap_CDIAC_data_2009.xls

% TODO: Does Marland et al. (2006) appear anywhere in this?
% Another option: give compact variable names, but be clear in defining
% them at the beginning, and comment a lot
% 9/16: last thing I did was line 38 (create new time vector)

function [ff1] = LoadFossilFuelData(timeStepPerYear)

%% load fossil fuel data

%Read in the data to be interpolated
ff=xlsread('BP_extrap_CDIAC_data_2009.xls'); %fossil fuel data matrix
ff_yr = ff(:,1); %year vector for fossil fuel data
ff_emis = ff(:,2); %emission vector for fossil fuel data

%Get rid of any NaN values
ff_emis(isnan(ff_yr)) = [];
ff_yr(isnan(ff_yr)) = [];

%Compute cumulative flux
cumFlux = 0; %cumulative flux counter
ff_emis_cum = ff_emis; %initialize cumulative fossil fuel vector
yr_cum = ff_yr; %initialize year vector for the eventual cumulative fossil fuel matrix

for i = 1:length(ff_yr) %Loop through year matrix to calculate the cumulative FF emissions for each year
    cumFlux = cumFlux + ff_emis(i);
    ff_emis_cum(i) = cumFlux; %set element in cumulative emission array to computer cumulative flux
    yr_cum(i) = ff_yr(i)+1;  %add 1 because integral is valid at the end of the calendar year
end

%Create new time array for interpolated data
yr_interp = ff_yr(1):(1/timeStepPerYear):ff_yr(end)+1;

%Do interpolation, creating matrix ff_interp to hold interpolated data
ff_interp = interp1(yr_cum,ff_emis_cum,yr_interp,'spline');

%TODO: what is happening here?
%Computed fluxes by taking discrete time derivative of integrated fluxes
fos = ff_interp;  %creates dummy arrays
yr = yr_interp;
%fos(length(yr_interp)) = []; %setting last element to empty?
%yr(length(yr_interp)) = []; %setting last element to empty?
for i = 1:(length(fos)-1) %loop through 
    fos(i) = timeStepPerYear.*(ff_interp(i+1)-ff_interp(i));
    yr(i) = yr_interp(i);
end

% convert to ppm

fosppm = fos/2.12;

% TODO: rename these output matrices-- differnentiate from fos, ff_interp 
% what are the differences between all of these matrices? Time deriv
% integrated fluxes etc
ff1(:,2) = fosppm(1,:); %what ils this doing?
ff1(:,1) = yr(1,:);
