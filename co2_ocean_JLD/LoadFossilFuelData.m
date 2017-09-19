% file LoadFossilFuelData.m
% formerly "load_fossil2.m"
%
% author Lauren Rafelski, modified by Julia Dohner
%
% brief Loads fossil fuel emissions data from BP Statistical Review of 
% World Energy (BP, 2008)
%
% Changes by LR:
% 3/23/09: Change input file to BP_extrap_CDIAC_data_2007.xls
% 1/10/11: Change input file to BP_extrap_CDIAC_data_2009.xls
%
% TODO: Does Marland et al. (2006) appear anywhere in this?
% Another option: give compact variable names, but be clear in defining
% them at the beginning, and comment a lot

function [fossilFuelData] = LoadFossilFuelData(timeStepPerYear)

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
yr_interp_intermed = ff_yr(1):(1/timeStepPerYear):ff_yr(end)+1;

%Do interpolation, creating matrix ff_interp to hold interpolated data
ff_interp_intermed = interp1(yr_cum,ff_emis_cum,yr_interp_intermed,'spline');

%Computed fluxes by taking discrete time derivative of integrated fluxes
ff_interp_flux = ff_interp_intermed;  %creates dummy arrays
for i = 1:(length(ff_interp_flux)-1) %loop through to fill ff_inter_flux 
    ff_interp_flux(i) = timeStepPerYear.*(ff_interp_intermed(i+1)-ff_interp_intermed(i));
end

% convert to ppm
fosppm = ff_interp_flux/2.12;

%Transpose year and emission values into columns 1 and 2 of output matrix,
%respectively
fossilFuelData(:,1) = yr_interp_intermed(1,:);
fossilFuelData(:,2) = fosppm(1,:);