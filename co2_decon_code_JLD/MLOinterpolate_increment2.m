% file MLOinterpolate_increment2.m
%
% author Lauren Rafelski, modified by Julia Dohner
% 
% inputs dtdelpCO2a,dpCO2a,year,dt,MLOSPOice_interp
% outputs [dtdelpCO2a,dpCO2a,year,dt,MLOSPOice_interp]
% 
% brief MLOinterpolate_increment2.m takes in a timestep and time 
% boundaries, uses merged CO2 data from Mauna Loa, South Pole Observatory 
% and Meure ice core records to return outputs for co2 increments for a 
% given timestep, rate of change of CO2, a year matrix, dt, and the merged 
% interpolated CO2 data. 
%
% See Figure 1.(A) of Rafelski (2008) for result of interpolation, merging of
% datasets.

%% 3/13/08: changed CO2 dataset to spline fit with sigma =0.6 to capture
%% 1940s plateau
%% 4/2/08: changed CO2 dataset to updated ice core data, sigma=0.6
%% 4/24/08: add in an option to decrease ice core data by 2 ppm; change to
%% linear interpolation for this option
%% 1/10/11: add updated CO2 dataset through 2010


%TODO: -ask Lauren what mlospo_meure, MLOSPOiceinterp, annincMLOSPO mean
% -make the for loop below way clearer (CO2 increment with monthly
% resolution, in ppm/year)
% -make all code fit onto 80 lines (linter for MATLAB?)

function [dtdelpCO2a,dpCO2a,year,dt,MLOSPOice_interp] = MLOinterpolate_increment2(timeStepPerYear,start_year,end_year)

%% load mlospo_meure matrix from co2_2011_2.mat file
%
% mlospo_meure is merged data from Mauna Loa, the South Pole Observatory,
% and ice core data before 1958 from MacFarling Meure et al. (2006)
%
load('co2_2011_2.mat', 'mlospo_meure'); 

dt = 1/timeStepPerYear; 

years0 = mlospo_meure(:,1); %fill years0 array with first column of mlospo_meure
CO2 = mlospo_meure(:,2); %fill CO2 array with second column of mlospo_meure

%% Create new time array

%copying year values (column 1 of mlospo_meure) by increment dt into years
%matrix

years = mlospo_meure(:,1):dt:mlospo_meure(end,1);

%% Do interpolation

%interpolate the CO2 data
CO2_interp = interp1(years0,CO2,years,'spline'); %spline is an option for data interpolation

% Merged and interpolated data (Mauna Loa, South Pole, Meure ice record)
MLOSPOice_interp(:,1) = years; % fill first column of MLOSPOiceinterp with years
MLOSPOice_interp(:,2) = CO2_interp; % fill second column with interpolated CO2 data

%% CO2 increment with monthly resolution, in ppm/year
% n = 7 is 7/1958, last value is 7/2005

%TODO: make this for loop way clearer
% I need A LOT of help on this section. Moving on for now.

%allocating space for annincMLOSPO matrix before for loop
%TODO: if timestep = 12, annincMLOSPO is 4430 long
dtdelpCO2a_numRows = length(((timeStepPerYear/2)+1):(length(years)-(timeStepPerYear/2)));
dtdelpCO2a = NaN(dtdelpCO2a_numRows,2);


%TODO: (ts=12, length(years) = 4442 -> "n = 7:(4442-6)")
for n = ((timeStepPerYear/2)+1):(length(years)-(timeStepPerYear/2))
    %fill first column of annincMLOSPO with first column of MLOSPOiceinterp
    dtdelpCO2a(n,1) = MLOSPOice_interp(n,1);
    %fill second column of annincMLOSPO with something complicated from
    %MLOSPOiceinter
    dtdelpCO2a(n,2) = MLOSPOice_interp(n+(timeStepPerYear/2),2) - MLOSPOice_interp(n-(timeStepPerYear/2),2);
end

year = start_year:dt:end_year;

%% Calculate change in atmospheric concentration

i = find(floor(100*MLOSPOice_interp(:,1)) == floor(100*(start_year+(1/24))));
dpCO2a(:,1) = MLOSPOice_interp(i:length(MLOSPOice_interp),1); 
dpCO2a(:,2) = MLOSPOice_interp(i:length(MLOSPOice_interp),2)-MLOSPOice_interp(i,2);
