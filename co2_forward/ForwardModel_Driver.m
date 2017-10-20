% file GlobalCO2_Forward_Driver.m
% 
% author Julia Dohner, borrowing from Lauren Rafelski
% 
% note Be sure to run defaults.m in outer folder before running this code
% 
% brief Fitting co2 record forward by 50 years, updating dpCO2a
% 
%
% subroutines: 
% HILDA response model r (calculates HILDA ocean uptake - A.2.2 in Joos 96)
% MLPulseResponse (calculates mixed layer uptake - Eq. 3, 6b in Joos 1996)
% biobox (Calculates co2 uptake of fast land box model (driven by co2
% fertilization))

% Oct 19, 2017 - trying to figure out what years the dpCO2a record
% currently spans. So it looks like the year vectors are useless in the
% original code. The axis that the lines are plotted on is just predefined
% and forced. There doesn't seem to be information in the year vector or
% the time column of the dpco2 vectors etc (i.e. they all seem to be from
% year 0005). But the range on the axis is 1800 to 2010. Super frustrating
% because I don't know which dpco2 values correspond to which dates. Agh
% just realized all of the above text was from only looking at the ocean
% code. Oh I see now. The time vectors are just years, not matlab
% timestamps. Assuming years are incremented by the increments in co2
% measurements

clear all

% give access to data files in co2_forward_data folder
addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));

LU = 1; %1 = high land use scenario; 2 = low land use scenario

nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?

filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered


load land_temp.mat % land temperature records

load npp_T.mat % NPP-weighted temperature record

load landwt_T_2011.mat % land temperature anomaly

% load CO2 sources and sinks
% landusemo: land use emissions, in ppm/year, interpolated to monthly
% resolution
% ff1: fossil fuel emissions
% fas: ocean flux per m^2
% Aoc: ocean surface area
% extratrop_landmo: extratropical land use emissions

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3; 

clear year start_year end_year ts

ts = 12; % timesteps per year
start_year = 1850;
end_year = 2009+(7/12); % is this now 2010?

beta = [0.5;2]; % initial guesses for model fit

% TODO: do we need the co2 record for the forward fitting? Yes
[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year); 
 
% % Extend land use record by making recent emissions equal to last
% % record
% landusemo(1874:1916,1) = year(1874:1916);
% landusemo(1874:1916,2) = landusemo(1873,2);
% % % 
% %Extend extratropical emissions by assuming emissions are zero
% extratrop_landmo(1802:1916,1) = landusemo(1802:1916,1);
% extratrop_landmo(1802:1916,2) = 0;
% 
% %% Calculate residual land uptake - DECONVOLUTION (calculates B) (land code)
% % run to 8/2009
% % using high land use emissions
% residual(:,1) = year(1,1:1916);
% residual(:,2) = dtdelpCO2a(2521:4436,2) - ff1(1189:3104,2)....
% + Aoc*fas(601:2516,2) - landusemo(1:1916,2);
% 
% % using extratropical emissions only
% residual2(:,1) = year(1,1:1916);
% residual2(:,2) = dtdelpCO2a(2521:4436,2) - ff1(1189:3104,2)....
% + Aoc*fas(601:2516,2) - extratrop_landmo(1:1916,2);


%function [t,r]= HILDAResponse(year)
%function [fas,dpCO2s]= MLPulseResponse(year,dpCO2a,c,h,kg,T,Aoc,r,dt)
%function [C1dt,C2dt,delCdt,delC1,delC2] = BioBoxResponse(eps,Q1a,Q2a,ts,year,dpCO2a,T)
%function [landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3;

%% THE MOTHERLOOP

end_year_forward = end_year + 50; % how far forward we'll be projecting


% I might as well loop through all the past data to see if I can fit it
% well, as well as looping through future times

% issue: dpco2a, fas, etc etc all seem to be slightly different lengths
% I suppose I don't have to initialize their lengths- I can just make
% vectors of the forward fits and fill those.... ugh not sure what to do
% here. I'm running the LR land code to find the size of all the matrices

% fas has two columns: year, and fas. year ends at 2010. Comes from
% getsourcesink

% dpCO2a 

%[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3; 

for i = start_year:dt:end_year
    %function [t,r]= HILDAResponse(year)
    %function [fas,dpCO2s]= MLPulseResponse(year,dpCO2a,c,h,kg,T,Aoc,r,dt)
    %function [C1dt,C2dt,delCdt,delC1,delC2] = BioBoxResponse(eps,Q1a,Q2a,ts,year,dpCO2a,T)
    %function [landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3;
end 
    

%% my own calculations/fiddling

% calculate timestep increment in LR code

diffArray = [];
for i = 1:(length(year)-1);
    diffArray(i) = year(i+1) - year(i);
end

% are all elements of diffarray = [0.0833333333332575]?

for i = 1:(length(year)-1);
    diffArray_test(i) = diffArray(i)-0.0833333333332575;
end

% conclusion: the timesteps are not *exactly* the same (they differ on the
% order of 10^(-16), but I think for our purposes we can say that the
% timesteps are consistent throughout. No! Be rigorous! Ugh. Cut off last 4
% decimal numbers? How are the years even calculated anyways?

% in MLOinterpolate_increment2: year = start_year:dt:end_year;

% THE TIMESTEP IS 1/DT WHICH IS JUST 1/12 = 0.0833 you're done there. The
% rest is weird roundoff errors


