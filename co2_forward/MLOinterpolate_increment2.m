% 3/13/08: changed CO2 dataset to spline fit with sigma =0.6 to capture
% 1940s plateau
% 4/2/08: changed CO2 dataset to updated ice core data, sigma=0.6
% 4/24/08: add in an option to decrease ice core data by 2 ppm; change to
% linear interpolation for this option
% 1/10/11: add updated CO2 dataset through 2010

function [annincMLOSPO,dpCO2a] = MLOinterpolate_increment2(ts,start_year,end_year)

dt = 1/ts;
    
%% processing data for co2 up to 2011 (in keeping with LR original code)
% loads in mlospo_meure (4427x2)
% starts 1640
% ends Feb 2010
% roughly monthly increments, but not perfectly. interpolate later
load co2_2011_2.mat

% years0 = mlospo_meure(:,1); % kind of unnecessary to copy over
% CO2_2011 = mlospo_meure(:,2);

% Create time array for co2_2011
year_2011 = mlospo_meure(1,1):1/ts:mlospo_meure(end,1);

% Do interpolation
% force to be exactly monthly resolution
co2_2011 = interp1(mlospo_meure(:,1),mlospo_meure(:,2),year_2011,'spline');

MLOSPOiceinterp_2011(:,1) = year_2011; % at monthly resolution
MLOSPOiceinterp_2011(:,2) = co2_2011;

%% everything past where co2_2011 ends

year_full = start_year:1/ts:end_year; 
% starts at year 1 (?) incremented by year
CO2_2016 = csvread('mergedCO2_2016.csv');
CO2_2016mo(:,1) = year_full;
CO2_2016mo(:,2) = (interp1(CO2_2016(:,1),CO2_2016(:,2),year_full)).';

% shorten to past end of co2_2011
i = find(floor(100*MLOSPOiceinterp_2011(:,1)) == floor(100*(start_year+(1/24))));
j = find(floor(100*(CO2_2016mo(:,1)+(1/24))) == floor(100*MLOSPOiceinterp_2011(end,1)));

co2_combine(:,1) = year_full; 
co2_combine(:,2) = [MLOSPOiceinterp_2011(i:end,2); CO2_2016mo(j+1:end,2)];


%% Calculate CO2 increment with monthly resolution, in ppm/year
% n = 7 is 7/1958, last value is 7/2005

for n = ((ts/2)+1):(length(co2_combine)-(ts/2))
    annincMLOSPO(n,1) = co2_combine(n,1);
    annincMLOSPO(n,2) = co2_combine(n+(ts/2),2) - co2_combine(n-(ts/2),2);
end

year = start_year:dt:end_year;

%% Calculate change in atmospheric concentration
% dpcCO2a = change in atmospheric co2 from preindustrial
%j = find(floor(100*co2_combine(:,1)) == floor(100*(start_year+(1/24))));
dpCO2a(:,1) = co2_combine(:,1); 
dpCO2a(:,2) = co2_combine(:,2)-co2_combine(1,2);


% %% 3/13/08: changed CO2 dataset to spline fit with sigma =0.6 to capture
% %% 1940s plateau
% %% 4/2/08: changed CO2 dataset to updated ice core data, sigma=0.6
% %% 4/24/08: add in an option to decrease ice core data by 2 ppm; change to
% %% linear interpolation for this option
% %% 1/10/11: add updated CO2 dataset through 2010
% 
% function [annincMLOSPO,dpCO2a,year,dt,MLOSPOiceinterp] = MLOinterpolate_increment2(ts,start_year,end_year)
% 
% load co2_2011_2.mat
% 
% dt = 1/ts; 
% 
% years0 = mlospo_meure(:,1);
% 
% CO2 = mlospo_meure(:,2);
% 
% %% Create new time array
% years = [years0(1):1/ts:years0(length(years0))];
% 
% %% Do interpolation
% CO22 = interp1(years0,CO2,years,'spline');
% 
% MLOSPOiceinterp(:,1) = years;
% MLOSPOiceinterp(:,2) = CO22;
% 
% %% CO2 increment with monthly resolution, in ppm/year
% % n = 7 is 7/1958, last value is 7/2005
% 
% for n = ((ts/2)+1):(length(years)-(ts/2))
%     annincMLOSPO(n,1) = MLOSPOiceinterp(n,1);
%     annincMLOSPO(n,2) = MLOSPOiceinterp(n+(ts/2),2) - MLOSPOiceinterp(n-(ts/2),2);
% end
% 
% year = start_year:dt:end_year;
% 
% %% Calculate change in atmospheric concentration
% i = find(floor(100*MLOSPOiceinterp(:,1)) == floor(100*(start_year+(1/24))));
% dpCO2a(:,1) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),1); 
% dpCO2a(:,2) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),2)-MLOSPOiceinterp(i,2);