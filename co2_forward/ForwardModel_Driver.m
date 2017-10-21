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

%end_year_forward = end_year + 50; % how far forward we'll be projecting


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

[t,r]= HILDAResponse(year); % r is kernel/pulse response fx

for i = start_year:dt:end_year
    
    %function [fas,dpCO2s]= MLPulseResponse(year,dpCO2a,c,h,kg,T,Aoc,r,dt)
    %function [C1dt,C2dt,delCdt,delC1,delC2] = BioBoxResponse(eps,Q1a,Q2a,ts,year,dpCO2a,T)
    %function [landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3;
    
    % ocean uptake
     %Calculate flux 
   fas(m,1) = year(m);
   fas(m,2) = (kg/Aoc)*(dpCO2a(m,2) - dpCO2s(m,2)); % air-sea flux of CO2
   
    w = conv(fas(1:m,2),r(1:m,2)); % convolve the air-sea flux and the pulse response function, as in Joos 1996

    % Calculate delDIC 
    delDIC(m+1,1) = year(m+1);
   
    delDIC(m+1,2) = (c/h)*w(m)*dt; % change in DIC

    %Calculate dpCO2s from DIC - from Joos 1996
    dpCO2s(m+1,2) = (1.5568 - (1.3993E-2)*T)*delDIC(m+1,2) + (7.4706-0.20207*T)*10^(-3)*...
        (delDIC(m+1,2))^2 - (1.2748-0.12015*T)*10^(-5)*(delDIC(m+1,2))^3 + (2.4491-0.12639*T)...
        *10^(-7)*(delDIC(m+1,2))^4 - (1.5468-0.15326*T)*10^(-10)*(delDIC(m+1,2))^5;
    
    
    
    % land uptake (biobox10)
    
      % fast box
    
      %dpco2a = 
      
    C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2)) - K1a*Q1a^((T(i,2)-T0)/10)*(C1 + delC1(i,2)); % temperature-dependent respiration
        
   % C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2))*(1 + Q1a*(T(i,2)-T0)) - K1a*(C1 + delC1(i,2)); % temperature-dependent photosynthesis
    
    % slow box
    
    C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2)) - K2a*Q2a^((T(i,2)-T0)/10)*(C2 + delC2(i,2)); % temperature-dependent respiration
        
   % C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2))*(1 + Q2a*(T(i,2)-T0)) - K2a*(C2 + delC2(i,2)); % temperature dependent photosynthesis  
    
    % box total change in concentrations
    
    delC1(i+1,2) = sum(C1dt(:,2))*dt;
    delC2(i+1,2) = sum(C2dt(:,2))*dt;
    
    % total flux into land
    
    delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
    
    B = delCdt;
    
    pco2a = pco2a + FF + LU - O - B; % updating pco2a
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


