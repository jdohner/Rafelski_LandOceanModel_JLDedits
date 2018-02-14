% forward model re-engineered

% file ForwardModel_Driver2.m
% version where data starts at 1850, ends at 2006 across all files
% re-engineered so that ocean and land uptake are calculated from observed
% atmospheric co2 record, but then sum the various fluxes (FF, LU, ocean
% and land) at the end to come up with new modeled atmospheric co2 record


clear all; %close all;

%% define run types

beta = [0.79;4.91]; % initial guesses for model fit
% VHM-V: [0.83;3.3]
% CHM-V: [0.79;4.91]

% predict = 1 --> prognostic, calculating dpCO2a in motherloop
% predict = 0 --> diagnostic (single deconvolution), feeding dpCO2a from
% MLOinterp, calculating residual land uptake in motherloop
predict = 1;

% define what kind of run you want to do
LU = 1; %1 = high land use scenario; 2 = low land use scenario
nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?
filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered

% nitrogen fertilization case, set epsilon = 0 to not allow CO2 fertilization
if nitrogen == 1 
    eps = 0;
    gamma = beta(1); % gamma is an input only for bioboxtwo_subN (N fertlz)
    Q1a = beta(2);
    Q2a = 1;
    
else % For CO2 fertilization model
    eps = beta(1);
    Q1a = beta(2);
    Q2a = 1;
end

%% define time frame 

start_year = 1800; % start year from LR ocean model
if predict == 1
    end_year = 2016;%2009+(7/12);%2016; %
else
    end_year = 2009+(7/12); % end year from LR land model
end

ts = 12; % number of data points/year, 12 in both land and ocean
dt = 1/ts;
year = start_year:dt:end_year;
year2 = year';


%% loading data

% give access to data files in co2_forward_data folder
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data'));
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data_2016'));



% loading temp data
addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));

load land_temp.mat % land temperature records
load npp_T.mat % NPP-weighted temperature record
load landwt_T_2011.mat % land temperature anomaly



%% constants

% ocean constants
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996


% land constants - define box sizes in ppm
Catm = 600/2.12; % around 283 ppm (preindustrial)
C1 = 110/2.12; % fast biosphere box, from old model
C2 = 1477/2.12; % slow biosphere box, changed to be same as 1 box

% Rate constants
 K1a = 1/2.5; 
 Ka1 = K1a*C1/Catm;

K2a = 1/60; % slow box to atmosphere
Ka2 = K2a*C2/Catm;

% To make temperature-independent: set Q1 and Q2 to 1
% TODO: address line above

%% initialize vectors

% dpCO2a is the change in atmospheric CO2 from preindustrial value
% dpCO2a = zeros(length(year),2); 
% dpCO2a(:,1) = year; 
dpCO2s = zeros(length(year),2); % dissolved CO2
dpCO2s(:,1) = year(:,1);
integral = zeros(length(year),length(year));
delDIC = zeros(length(year),2); 
delC1(:,1) = year2(:,1);
delC1(:,2) = zeros(size(year2));
delC2(:,1) = year2(:,1);
delC2(:,2) = zeros(size(year2));
delCdt(:,1) = year2(:,1);
C1dt(:,1) = year2(:,1);
C2dt(:,1) = year2(:,1);
delC1(length(year2)+1,1) = year2(length(year2),1)+dt;
delC2(length(year2)+1,1) = year2(length(year2),1)+dt;
residualLandUptake = [];
residualLandUptake(:,1) = year2;
% dtdelpCO2a = []; 
% dtdelpCO2a(:,1) = year2;
fas = zeros(length(year2),2);
fas(:,1) = year2;
integrationCheck(:,1) = year2;
integrationCheck(:,2) = 0;
cum_ff(:,1) = year2;
cum_ff(:,2) = 0;
cum_lu(:,1) = year2;
cum_lu(:,2) = 0;
cum_ocean(:,1) = year2;
cum_ocean(:,2) = 0;
cum_land(:,1) = year2;
cum_land(:,2) = 0;
dpCO2a_calc(:,1) = year2;
dpCO2a_calc(:,2) = 0;
dtdelpCO2a_calc(:,1) = year2;
dtdelpCO2a_calc(:,2) = 0;


%% set up ocean

% Response function to calculate ocean uptake
[t,r] = HILDAResponse(year);


%% set up land

[ff, landuse, landuseExtra] = getSourceSink2(predict,year);

% tland4: begins in 1880. tland4 extends record back to 1800 using the mean
% temperature for the first years
% what's the difference between tland4 and landtglob?


[temp_anom, T0] = tempRecord(tland4, landtglob, dt, year, end_year);



%% Get observed CO2 record

% MLOinterp isn't called anywhere else in the code
[dtdelpCO2a_obs,dpCO2a_obs,year_obs,dt_obs,CO2a_obs] = ... 
    MLOinterpolate_increment2(ts,start_year,end_year); 

% add on recent co2 data
CO2_2016 = csvread('mergedCO2_2016.csv');
CO2_2016mo(:,1) = year;
CO2_2016mo(:,2) = (interp1(CO2_2016(:,1),CO2_2016(:,2),year)).';
co2_preind = mean(CO2_2016(1:1000,2));
CO2_2016mo(:,2) = CO2_2016mo(:,2) - co2_preind;
i = find(CO2_2016mo(:,1) == 2010); % come back to this with actual indexing

newCO2obs(:,1) = year2;
newCO2obs(:,2) = 0;
newCO2obs(1:2520,2) = dpCO2a_obs(1:2520,2);
newCO2obs(i:end,2) = CO2_2016mo(i:end,2);
dpCO2a_obs = newCO2obs;


% co2_predict = dpCO2a;
% co2_predict(:,2) = co2_predict(:,2) + co2_preind;

% this is bootleg, but it'll have to do for now:
% ALWAYS BE WARY OF THIS WHEN DEBUGGING:
dtdelpCO2a_obs = dtdelpCO2a_obs(1919:4434,:); % indices of closest values 
% to 1800 and 2006


% not proper but it'll have to do for now:
%dpCO2a_obs = dpCO2a_obs(1:2516,:);

% cut off dpco2a to end at end_year
%     co2_end = find(dpCO2a_obs(:,1) == year(end)); % buggy line - need data thru 2016
%     dpCO2a_obs = dpCO2a_obs(1:co2_end,:);
%     

%[c index] = min(abs(N-V(1)))

% this doesn't quite work since closest value is below (floor problem)
% i = find(floor(100*dtdelpCO2a_obs(:,1)) == floor(100*start_year));
% j = find(floor(100*dtdelpCO2a_obs(:,1)) == floor(100*end_year));
% dtdelpCO2a_obs = dtdelpCO2a_obs(i:j,:);

%% Motherloop v2.0

for i = 1:length(year2) 
    
    fas(i,1) = year(i);
    fas(i,2) = (kg/Aoc)*(dpCO2a_obs(i,2) - dpCO2s(i,2)); % air-sea flux of CO2

    % convolve the air-sea flux and the pulse response function (Joos 1996)
    w = conv(fas(1:i,2),r(1:i,2)); 

    if i < length(year2)
    % Calculate delDIC
        delDIC(i+1,1) = year(i+1); 
        delDIC(i+1,2) = (c/h)*w(i)*dt; % change in DIC

    %Calculate dpCO2s from DIC - from Joos 1996
    dpCO2s(i+1,2) = (1.5568 - (1.3993E-2)*T_const)*delDIC(i+1,2) + ...
        (7.4706-0.20207*T_const)*10^(-3)*(delDIC(i+1,2))^2 - ...
        (1.2748-0.12015*T_const)*10^(-5)*(delDIC(i+1,2))^3 + ...
        (2.4491-0.12639*T_const)*10^(-7)*(delDIC(i+1,2))^4 - ...
        (1.5468-0.15326*T_const)*10^(-10)*(delDIC(i+1,2))^5;
    end
    
    % land uptake - code from bioboxtwo_sub10_annotate.m for CO2
    
    % fertilization, from bioboxtwo_subN_annotate.m for nitrogen fert
    
    if nitrogen == 1 % nitrogen fertilization case
        
        % fast box
        % temperature-dependent respiration
        C1dt(i,2) = Ka1*Catm*(1 + eps*dpCO2a_obs(i,2)/Catm + gamma*ff(i+a-1,2))...
            - K1a*Q1a^((temp_anom(i,2)-T0)/10)*(C1 + delC1(i,2)); 
        % temperature-dependent photosynthesis
        % C1dt(m,2) = Ka1*Catm*(1 + eps*dpCO2a_obs(m,2)/Catm + gamma*ff1(m+a-1,2))...
        %*(1 + Q1a*(temp_anom(m,2)-T0)) - K1a*(C1 + delC1(m,2)); 
        
        % slow box
        % temperature-dependent respiration 
        C2dt(i,2) = Ka2*Catm*(1 + eps*dpCO2a_obs(i,2)/Catm + gamma*ff(i+a-1,2))...
            - K2a*Q2a^((temp_anom(i,2)-T0)/10)*(C2 + delC2(i,2));
        % temperature-dependent photosynthesis
        % C2dt(m,2) = Ka2*Catm*(1 + eps*dpCO2a_obs(m,2)/Catm + gamma*ff1(m+a-1,2))...
        %*(1 + Q2a*(temp_anom(m,2)-T0)) - K2a*(C2 + delC2(m,2)); 
        
        if i < length(year2)
            % box total change in concentrations
            delC1(i+1,2) = sum(C1dt(:,2))*dt;
            delC2(i+1,2) = sum(C2dt(:,2))*dt;
        end

        % total flux into land
        delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
    
    else % co2 fertilization case (nitrogen == 0)

        % land uptake (biobox10)
        % fast box
        % temperature-dependent respiration 
        C1dt(i,2) = Ka1*(Catm + eps*dpCO2a_obs(i,2)) - K1a*Q1a^((temp_anom(i,2)-T0)/10)...
            *(C1 + delC1(i,2)); 
        % temperature-dependent photosynthesis
        % C1dt(i,2) = Ka1*(Catm + eps*dpCO2a_obs(i,2))*(1 + Q1a*(temp_anom(i,2)-T0))...
        %- K1a*(C1 + delC1(i,2)); 

        % slow box
        % temperature-dependent respiration 
        C2dt(i,2) = Ka2*(Catm + eps*dpCO2a_obs(i,2)) - K2a*Q2a^((temp_anom(i,2)-T0)/10)...
            *(C2 + delC2(i,2)); 
        % temperature dependent photosynthesis  
        % C2dt(i,2) = Ka2*(Catm + eps*dpCO2a_obs(i,2))*(1 + Q2a*(temp_anom(i,2)-T0))...
        % - K2a*(C2 + delC2(i,2)); 
        
        if i < length(year2)
            % box total change in concentrations
            delC1(i+1,2) = sum(C1dt(:,2))*dt;
            delC2(i+1,2) = sum(C2dt(:,2))*dt;
        end
        
        % total flux into land
        delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
        
    end
    
    % calculated modeled atmospheric co2
    
    dtdelpCO2a_calc(i,2) =  ff(i,2) - Aoc*fas(i,2) - delCdt(i,2) + landuse(i,2);
        
    if i < length(year2)
            dpCO2a_calc(i+1,2) = dpCO2a_calc(i,2) + dtdelpCO2a_calc(i,2)/12; 
    end
    
    %         % For debugging:
%         %     fas(i,2) = 0;
%              dtdelpCO2a(i,2) = 0;
%              ff1(i,2) = 0;
%              landusemo(i,2) = 0;
%              B(i,2) = 0;

        cum_ff(i,2) = sum(ff(1:i,2)/12);
        cum_lu(i,2) = sum(landuse(1:i,2)/12);
        cum_ocean(i,2) = sum(-Aoc*fas(1:i,2)/12); % negative sign here because JLD convention, 
        cum_land(i,2) = sum(-delCdt(1:i,2)/12);  % neg sign here bc JLD not yet applied


        integrationCheck(i,2) =  cum_ff(i,2)  + cum_ocean(i,2) + cum_land(i,2) + cum_lu(i,2);

end

% update sign convention
fas(:,2) = -1*fas(:,2);
delCdt(:,2) = -1*delCdt(:,2);

%% plotting

ss = get(0,'screensize');
width = ss(3);
height = ss(4);


    
    % annual mass balance check
    annualMB = [];
    annualMB(:,1) = year2; 
    annualMB(:,2) = ff(:,2)  + Aoc*fas(:,2) + delCdt(:,2)...
            - dtdelpCO2a_calc(:,2) + landuse(:,2);
    
    figure('name','annual mass balance - PREDICT');
    plot(annualMB(:,1),annualMB(:,2));
    axis([1800 2010 -10 10])
    title('mass balance = atmos - land - ff + Aoc*fas')
    xlabel('Year')
    
    % sources and sinks (to compare to LR)
    x = 0*year;
    H = figure('name','sources and sinks - PREDICT');
    plot(delCdt(:,1),delCdt(:,2),'-g',ff(:,1), ff(:,2), ...
        '-k', dtdelpCO2a_calc(:,1),dtdelpCO2a_calc(:,2),'-r', fas(:,1),Aoc*fas(:,2),'-b',year(1,:),x,'--k');
    axis([1800 2010 -10 10])
    legend('land flux','fossil fuel','modeled atmosphere','ocean','Location','SouthWest')
    title('sources and sinks - PREDICT')
    xlabel('Year ')
    ylabel('ppm/year  Positive = source, negative = sink ')


    % cumulative mass balance
    figure('name','Cumulative Mass Balance');
    diff = integrationCheck(:,2)-dpCO2a_obs(:,2);
    plot(dpCO2a_calc(:,1), dpCO2a_calc(:,2), ...
        ff(:,1), cum_ff(:,2), ff(:,1), cum_ocean(:,2), ff(:,1),cum_land(:,2),...
        ff(:,1),cum_lu(:,2),ff(:,1),diff(:,1),year(1,:),x,'--k');
    legend('dpCO2a_calc','Cumulative FF','Cumulative ocean sink', 'Cumulative land', 'cumulative landuse','Cumulative mass balance','Location','NorthWest')
    set(findall(gca, 'Type', 'Line'),'LineWidth',4);
    
    % compare predicted co2 to observed - basically fig 7
    
    CO2_2016 = csvread('mergedCO2_2016.csv');
    CO2_2016mo(:,1) = year;
    CO2_2016mo(:,2) = (interp1(CO2_2016(:,1),CO2_2016(:,2),year)).';
    figure('name','CO2: Predicted vs. Observed');
    co2_preind = mean(CO2_2016(1:1000,2));
    co2_predict = dpCO2a_calc;
    co2_predict(:,2) = co2_predict(:,2) + co2_preind;
    co2_diff = CO2_2016mo(:,2)-co2_predict(:,2);

    % from paper: add integration constant to achieve best agreement with
    % co2 record from 1959 to 1979
    meandiff = mean(co2_diff(1909:2149)); % mean difference over 1959-1979
    predict2 = co2_predict(:,2)+meandiff;
    %plot(CO2_2016mo(:,1), CO2_2016mo(:,2),co2_predict(:,1),co2_predict(:,2),CO2_2016mo(:,1),co2_diff);
    %  fixed fraction of fossil fuel emissions (fit to observations) -- fit
    %  between 1959 and 1979
    ff57(:,1) = cum_ff(:,1);
    ff57(:,2) = cum_ff(:,2)*0.57+co2_preind;
    co2_diff2 = CO2_2016mo(:,2)-ff57(:,2); % getting factor by which to move up 57% ff
    meandiff2 = mean(co2_diff2(1909:2149)); % mean difference over 1959-1979
    ff57_2 = ff57(:,2)+meandiff2;
    
    plot(CO2_2016mo(:,1), CO2_2016mo(:,2),co2_predict(:,1),predict2, ff57(:,1),ff57_2, '--k');
    legend('Observations','Temperature-dependent model (CHM-V)','Fixed fraction of ff emissions (fit to observations');
    set(findall(gca, 'Type', 'Line'),'LineWidth',4);
    xlim([1850 year(end)])
    ylim([275 400])
    
    figure('name','Atmospheric CO2 History - Figure 7, Rafelski (2009)');
    newdiff = predict2 - CO2_2016mo(:,2);
    plot(CO2_2016mo(:,1), newdiff, year(1,:),x,'--k');
    legend('predicted - observed co2');
    set(findall(gca, 'Type', 'Line'),'LineWidth',4);
    xlim([1850 year(end)])
    
    % recreating figure 7a
    figure('name','cumulative mass balance');
    diff = integrationCheck(:,2)-dpCO2a_obs(:,2);
    plot(dpCO2a_calc(:,1), dpCO2a_calc(:,2), ...
        ff(:,1), cum_ff(:,2), ff(:,1), cum_ocean(:,2), ff(:,1),cum_land(:,2),...
        ff(:,1),cum_lu(:,2),ff(:,1),diff(:,1),year(1,:),x,'--k');
    legend('dpCO2a_calc','Cumulative FF','Cumulative ocean sink', 'Cumulative land', 'cumulative landuse','Cumulative mass balance','Location','NorthWest')
    set(findall(gca, 'Type', 'Line'),'LineWidth',4);
    
    
    
% position figures
    vert = 300; %300 vertical pixels
    horz = 600; %600 horizontal pixels
    set(H,'Position',[(3*width/4)-horz/2, (height/2)-vert/2, horz, vert]);
    I = openfig('LR_plot.fig');
    set(I,'Position',[(width/4)-horz/2, (height/2)-vert/2, horz, vert]);


 