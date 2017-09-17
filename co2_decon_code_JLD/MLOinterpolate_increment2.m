% file MLOinterpolate_increment2.m
% author Lauren Rafelski, modified by Julia Dohner
% brief Defines function MLOinterpolate_increment2

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

function [annincMLOSPO,dpCO2a,year,dt,MLOSPOiceinterp] = MLOinterpolate_increment2(timeStepPerYear,start_year,end_year)

%load mlospo_meure matrix from co2_2011_2.mat file
load('co2_2011_2.mat', 'mlospo_meure'); 

dt = 1/timeStepPerYear; 

years0 = mlospo_meure(:,1); %fill years0 array with first column of mlospo_meure
CO2 = mlospo_meure(:,2); %fill CO2 array with second column of mlospo_meure

%% Create new time array
%copying year values (column 1 of mlospo_meure) by increment dt into years
years = mlospo_meure(:,1):dt:mlospo_meure(end,1);

%% Do interpolation
CO2_interp = interp1(years0,CO2,years,'spline'); %spline is an option for data interpolation

MLOSPOiceinterp(:,1) = years;
MLOSPOiceinterp(:,2) = CO2_interp;

%% CO2 increment with monthly resolution, in ppm/year
% n = 7 is 7/1958, last value is 7/2005

%TODO: make this for loop way clearer

%allocating space for annincMLOSPO matrix before for loop
annincMLOSPO_numRows = length(((timeStepPerYear/2)+1):(length(years)-(timeStepPerYear/2)));
annincMLOSPO = NaN(annincMLOSPO_numRows,2);

for n = ((timeStepPerYear/2)+1):(length(years)-(timeStepPerYear/2))
    %fill first column of annincMLOSPO with first column of MLOSPOiceinterp
    annincMLOSPO(n,1) = MLOSPOiceinterp(n,1); %annual increase Mauna Loa SPO?
    %fill second column of annincMLOSPO with something complicated from
    %MLOSPOiceinter
    annincMLOSPO(n,2) = MLOSPOiceinterp(n+(timeStepPerYear/2),2) - MLOSPOiceinterp(n-(timeStepPerYear/2),2);
end

year = start_year:dt:end_year;

%% Calculate change in atmospheric concentration
i = find(floor(100*MLOSPOiceinterp(:,1)) == floor(100*(start_year+(1/24))));
dpCO2a(:,1) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),1); 
dpCO2a(:,2) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),2)-MLOSPOiceinterp(i,2);
