% file ForwardModel_Driver2.m

clear all;

% give access to data files in co2_forward_data folder
addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));

%% set up ocean

% ~~~~~~~copying from jooshildascale_annotate2.m from here down~~~~~~~

ts = 12; % number of data points/year
start_year_ocean = 1800;
end_year_ocean = 2010;
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996

% don't get atmospheric co2 record from MLOinterp
[ff1] = load_fossil2(ts); % get fossil fuel emissions

% initializing year vector (from MLOinterpolate_increment2.m)
dt = 1/ts; 
year = start_year_ocean:dt:end_year_ocean;

% back to jooshildascale_annotate2.m

% Response function to calculate ocean uptake
[t,r] = HILDAResponse(year);

% in the original code, would call joos_general_fast_annotate2 here, but
% I'm putting part of that function into the motherloop, so here I'm just
% calling the prep lines from that function (pre-loop)

% dpCO2a is initialized in MLOinterp, but we need it for its dimensions to
% initialize dpCO2s
% dpCO2a as outputted from MLOinterp is a 2522x2 double
dpCO2a = [1:2522,2];
% these dimensions are set by the length of the record in mlospo_meure and
% by the timestep ts

dpCO2s = zeros(length(dpCO2a),2); % dissolved CO2
dpCO2s(:,1) = dpCO2a(:,1);
%fas = zeros(length(year),2); % removed because is output by get source
%sink
integral = zeros(length(year),length(year));
delDIC = zeros(length(year),2);

% the rest (lines 41 onward) in the ocean driver
% (jooshildascale_annotate2.m) just calculates residual land flux and plots

% TODO: make sure include lines 31-40 from joos_general_fast_annotate2.m
% after the motherloop

%% set up land

% ~~~~~~ copying from nonlin_land_Qs_annotate.m ~~~~~~~

% define what kind of run you want to do
LU = 1; %1 = high land use scenario; 2 = low land use scenario
nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?
filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered

% loading temp data
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

[landusemo,ff1,fas,extratrop_landmo] = getsourcesink_scale3; 

%ts = 12; % took out because already defined in ocean
start_year_land = 1850; % change to start_year_land
end_year_land = 2009+(7/12); % change to end_year_land

beta = [0.5;2]; % initial guesses for model fit

% leaving out a call to MLO interpolate to get dtdelpCO2a,dpCO2a,year,dt,CO2a
% we get a year vector from MLOinterp, but same as line 26 so won't repeat

% TODO: the following (extending records) can all happen in getsourcesink_scale3, but wait to
% change until can get the code running and working

% Extend land use record by making recent emissions equal to last
% record
landusemo(1874:1916,1) = year(1874:1916);
landusemo(1874:1916,2) = landusemo(1873,2);
% % 
%Extend extratropical emissions by assuming emissions are zero
extratrop_landmo(1802:1916,1) = landusemo(1802:1916,1);
extratrop_landmo(1802:1916,2) = 0;

% Note: when get to find model fit section, don't need to call nlinfit, this is
% just where do motherloop. Want to call bioboxes for their output, but
% don't need to do any fitting

% Question: Do I need to calculate the 10-year boxcar average of the residual
% land uptake? This whole next section:

% % calculate a 10-year running boxcar average of the residual land uptake
% % don't use 10 year mean before 1957, because data are already smoothed
% % (ice core)
% if(LU==1) %high land use
% [residual10] = l_boxcar(residual,10,12,1,length(residual),1,2);
% [residual10a] = l_boxcar(residual,10,12,1225,length(residual),1,2);
% residual10(1:1284,:) = residual(1:1284,:);
% residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);
% 
% decon = residual;
% elseif(LU ==2) % low land use
% [residual10] = l_boxcar(residual2,10,12,1,length(residual2),1,2);
% [residual10a] = l_boxcar(residual2,10,12,1225,length(residual2),1,2);
% residual10(1:1284,:) = residual2(1:1284,:);
% residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);
% 
% decon = residual2;
% end

% LR runs nlinfit to find model fit, then looks at covariances and
% correlations between model results and calculated land uptake, then gets
% uncertainties of best fit values. Then she redefines values of epsilon,
% gamma and Q1 based on the results. I'll instead use the values of
% epsilon, gamma and Q1 from Table 2 of her paper.

%% Redefine values of epsilon, gamma and Q1
betahat = zeros(2);

if(nitrogen == 1)
epsilon = 0;% betahat(2)
gamma = betahat(1);
Q1 = betahat(2);
Q2 = 1;%betahat(3)
else
epsilon = betahat(1); %0.79;%
Q1 = betahat(2); %4.91;
Q2 = 1; %betahat(2)
end

%% probably time for THE MOTHERLOOP

%before this call, year is a 1x2521 double vector (one row vector)
year2 = year';
% all along in this code, the year vector is 


% Question/TODO: Issue: the loop boundaries for ocean and land uptake
% (biobox, joos_general_fast_annotate2.m)
%
% at the moment, the ocean year vector is longer than the land year vector.
%
% for biobox, year vector is (defined in nonlin_land_Qs_annotate.m):
% ts = 12; % timesteps per year
% start_year = 1850;
% end_year = 2009+(7/12); 
% year = start_year:dt:end_year; (MLOinterp)
%
% for ocean (joos_general), year vector is:
% start_year = 1800;
% end_year = 2010;
% dt = 1/ts (MLOinterp)
% year = start_year:dt:end_year; (MLOinterp)

% note: fas goes from 1800 to 2010

% Index exceeds matrix dimensions.
% 
% Error in ForwardModel_Driver2 (line 173)
%     fas(i,2) = (kg/Aoc)*(dpCO2a(i,2) - dpCO2s(i,2)); % air-sea flux of CO2


% everything below is just the loop contents from
% bioboxtwo_sub10_annotate.m and joos_general_fast_annotate2, but lacks the
% stuff from both of those functions that help set up the loop.

% i loops through 2521 values

% biobox_sub10 loop setup:







for i = 1:length(year2);
    
    % ocean uptake - code from joos_general_fast_annotate2
    
    %Calculate flux 
    fas(i,1) = year(i);
    fas(i,2) = (kg/Aoc)*(dpCO2a(i,2) - dpCO2s(i,2)); % air-sea flux of CO2
   
    w = conv(fas(1:i,2),r(1:i,2)); % convolve the air-sea flux and the pulse response function, as in Joos 1996

    % Calculate delDIC 
    delDIC(i+1,1) = year(i+1); % filling time column for delDIC
    delDIC(i+1,2) = (c/h)*w(i)*dt; % change in DIC

    %Calculate dpCO2s from DIC - from Joos 1996
    dpCO2s(i+1,2) = (1.5568 - (1.3993E-2)*T_const)*delDIC(i+1,2) + (7.4706-0.20207*T_const)*10^(-3)*...
        (delDIC(i+1,2))^2 - (1.2748-0.12015*T_const)*10^(-5)*(delDIC(i+1,2))^3 + (2.4491-0.12639*T_const)...
        *10^(-7)*(delDIC(i+1,2))^4 - (1.5468-0.15326*T_const)*10^(-10)*(delDIC(i+1,2))^5;
    
    
    % land uptake - code from bioboxtwo_sub10_annotate.m for CO2
    % fertilization, from bioboxtwo_subN_annotate.m for nitrogen
    % fertilization
    
    % Question/TODO: what are we doing about temp-dependent photosynthesis
    % and respiration? TODO part is to make if cases for these
    if nitrogen == 1 % nitrogen fertilization case
        
        % fast box
        C1dt(i,2) = Ka1*Catm*(1 + eps*dpCO2a(i,2)/Catm + gamma*ff1(i+a-1,2)) - K1a*Q1a^((T(i,2)-T0)/10)*(C1 + delC1(i,2)); % temperature-dependent respiration
        % C1dt(m,2) = Ka1*Catm*(1 + eps*dpCO2a(m,2)/Catm + gamma*ff1(m+a-1,2))*(1 + Q1a*(T(m,2)-T0)) - K1a*(C1 + delC1(m,2)); % temperature-dependent photosynthesis

        % slow box
        C2dt(i,2) = Ka2*Catm*(1 + eps*dpCO2a(i,2)/Catm + gamma*ff1(i+a-1,2)) - K2a*Q2a^((T(i,2)-T0)/10)*(C2 + delC2(i,2)); % temperature-dependent respiration 
        % C2dt(m,2) = Ka2*Catm*(1 + eps*dpCO2a(m,2)/Catm + gamma*ff1(m+a-1,2))*(1 + Q2a*(T(m,2)-T0)) - K2a*(C2 + delC2(m,2)); % temperature-dependent photosynthesis

        % box total change in concentrations
        delC1(i+1,2) = sum(C1dt(:,2))*dt;
        delC2(i+1,2) = sum(C2dt(:,2))*dt;

        % total flux into land
        delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
    
    else % co2 fertilization case (nitrogen == 0)

        % land uptake (biobox10)
        % fast box
        % dpco2a = 
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
        
    end
    
    B = delCdt;
    
    pco2a = pco2a + FF + LU - O - B; % updating pco2a
end 

% post-loop code from joos_general_fast_annotate2.m
% calculate the flux for the last time point
fas(length(year),1) = year(length(year));
fas(length(year),2) = (kg/Aoc)*(dpCO2a(length(year),2) - dpCO2s(length(year),2));

zero1 = find(delDIC(:,1) == 0);
delDIC(zero1,2) = NaN;
delDIC(zero1,1) = NaN;

zero2 = find(dpCO2s(:,2) == 0);
dpCO2s(zero2,2)  = NaN;


%% Misc from nonlin_land_Qs_annotate.m

% Leaving all of this out. Question: probably don't need, right?
% % Run the best fit values (parameters from Table 2) in the model again to plot
% if nitrogen == 1
%     [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN_annotate(epsilon,Q1,Q2,gamma,ff1,ts,year2,dpCO2a,X);
% else 
%     [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 
% end
%     
% delCdt(:,2) = -delCdt(:,2);
% 
% % 10 year moving boxcar average of model result
% [delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);

%% Do "reverse deconvolution" to calculate modeled atmospheric change in
% CO2
% Question: Do I want to keep this?
% if(LU==1)
% newat(:,1) = year(1,1:1916);
% newat(:,2) =  ff1(1189:3104,2)....
% - Aoc*fas(601:2516,2) + landusemo(1:1916,2) + delCdt(:,2) ;
% 
% elseif(LU==2)
% newat(:,1) = year(1,1:1867);
% newat(:,2) =  ff1(1189:3055,2)....
% - Aoc*fas(601:2467,2) + extratrop_landmo(1:1867,2) + delCdt(:,2) ;
% end