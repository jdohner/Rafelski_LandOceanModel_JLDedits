% file ForwardModel_Driver2.m
% version where data starts at 1850, ends at 2006 across all files


clear all; %close all;

% give access to data files in co2_forward_data folder
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data_2016'));

% predict = 1 --> prognostic, calculating dpCO2a in motherloop
% predict = 0 --> diagnostic (single deconvolution), feeding dpCO2a from
% MLOinterp, calculating residual land uptake in motherloop
predict = 1;

% choose your sign convention
% JLD = 1: JLD sign convention (fas is positive into atmosphere)
% JLD = 0; LR sign convention (fas is positive into the ocean)
JLD = 1; 


%% set up ocean

start_year_ocean = 1800; % start year from LR ocean model
end_year_ocean = 2009+(7/12); % end year from LR land model
ts = 12; % number of data points/year, 12 in both land and ocean
dt = 1/ts;
year_ocean = start_year_ocean:dt:end_year_ocean;


% constants
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996

% Response function to calculate ocean uptake
[t,r] = HILDAResponse(year_ocean);

% dpCO2a is the change in atmospheric CO2 from preindustrial value
dpCO2a = zeros(length(year_ocean),2); 
dpCO2a(:,1) = year_ocean; 

dpCO2a_LRdata = load('dpCO2a_LRocean.mat');
dpCO2s_LRdata = load('dpCO2s_LRocean.mat');
delDIC_LRdata = load('delDIC_LRocean.mat');
w_LRdata = load('w_LRocean.mat');
dpCO2a_LRocean = dpCO2a_LRdata.dpCO2a;
dpCO2s_LRocean = dpCO2s_LRdata.dpCO2s;
delDIC_LRocean = delDIC_LRdata.delDIC;
w_LRocean = w_LRdata.w;

dpCO2a(1:10,2) = dpCO2a_LRocean(1:10,2);

dpCO2s = zeros(length(dpCO2a),2); % dissolved CO2
dpCO2s(:,1) = dpCO2a(:,1);
integral = zeros(length(year_ocean),length(year_ocean));
delDIC = zeros(length(year_ocean),2);

% TODO: make sure include lines 31-40 from joos_general_fast_annotate2.m
% after the motherloop

%% set up land

% ~~~~~~ copying from nonlin_land_Qs_annotate.m ~~~~~~~

% define what kind of run you want to do
LU = 1; %1 = high land use scenario; 2 = low land use scenario
nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?
filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered

% loading temp data
addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));

load land_temp.mat % land temperature records
load npp_T.mat % NPP-weighted temperature record
load landwt_T_2011.mat % land temperature anomaly

% create year vector for land to be called in landusemo calls
start_year_land = start_year_ocean; 
end_year_land = end_year_ocean; 
year_land = start_year_land:dt:end_year_land;

[ff1,landusemo,extratrop_landmo] = getsourcesink_scale4(predict,start_year_ocean,end_year_ocean,year_land); % ff1 load land



beta = [0.5;2]; % initial guesses for model fit


%%%%
% getsourcesink outputs end at 2006
% Extend land use record by making recent emissions equal to last
% record
j = length(landusemo);
k = length(year_ocean);
extendYears = year_ocean(1,j+1:k);
%landusemo(:,1) = [landusemo(:,1), extendYears];

landusemo(j+1:k,1) = year_ocean(1,j+1:k);
landusemo(j+1:k,2) = landusemo(j,2);

%Extend extratropical emissions by assuming emissions are zero
extratrop_landmo(j+1:k,1) = landusemo(j+1:k,1);
extratrop_landmo(j+1:k,2) = 0;


% TODO: the following (extending records) can all happen in 
% getsourcesink_scale3, but wait to
% change until can get the code running and working

%Extend extratropical emissions by assuming emissions are zero
% extratrop_landmo(1802:1916,1) = landusemo(1802:1916,1);
% extratrop_landmo(1802:1916,2) = 0;

% tland4: begins in 1880. tland4 extends record back to 1800 using the mean
% temperature for the first years

% do a moving boxcar average of the land temperature: 1 year average
[avg_temp] = l_boxcar(tland4,1,12,1,2483,1,2); 

avg_temp(1:6,2) = avg_temp(7,2); % make the first 6 points 
% Question: why does avg_temp has zeros for first 6 time points?
% do these time points just not get called? only the corresponding values?

temp_anom(1:606,1) = avg_temp(1:606,1); % fill time column of temp_anom
temp_anom(1:606,2) = landtglob(1,2); %landtglob are temp anomalies
temp_anom(607:2516,1) = landtglob(1:1910,1);
temp_anom(607:2516,2) = landtglob(1:1910,2);

clear avg_temp;

%% OBSERVATION DIAGNOSTIC DEBUGGING SECTION

% MLOinterp isn't called anywhere else in the code
[dtdelpCO2a_obs,dpCO2a_obs,year_obs,dt_obs,CO2a_obs] = ... 
    MLOinterpolate_increment2(ts,start_year_ocean,end_year_ocean); 

% truncate obs output to match my timeframe
%i = find(floor(100*MLOSPOiceinterp(:,1)) == floor(100*(start_year+(1/24))));
% startIndex = find(floor(dtdelpCO2a_obs(:,1) == start_year_ocean);
% endIndex = find(dtdelpCO2a_obs(:,1) == end_year_ocean);
% dtdelpCO2a_obs = dtdelpCO2a_obs(startIndex:endIndex,:);
% startIndex = 1919; % this is the closest date to 1800- it's actually 
% halfway through december of 1799

% this is bootleg, but it'll have to do for now:
% ALWAYS BE WARY OF THIS WHEN DEBUGGING:
dtdelpCO2a_obs = dtdelpCO2a_obs(1919:4434,:); % indices of closest values 
% to 1800 and 2006


%% probably time for THE MOTHERLOOP

%before this call, year is a 1x2521 double vector (one row vector)
year_ocean2 = year_ocean';

% biobox_sub10 loop setup:

T0 = temp_anom(1,2);

% Define box sizes in ppm

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



year_land_trans = year_land';
delC1(:,1) = year_land_trans(:,1);
delC1(:,2) = zeros(size(year_land_trans));
delC2(:,1) = year_land_trans(:,1);
delC2(:,2) = zeros(size(year_land_trans));
delCdt(:,1) = year_land_trans(:,1);
C1dt(:,1) = year_land_trans(:,1);
C2dt(:,1) = year_land_trans(:,1);
delC1(length(year_land_trans)+1,1) = year_land_trans(length(year_land_trans),1)+dt;
delC2(length(year_land_trans)+1,1) = year_land_trans(length(year_land_trans),1)+dt;



% this is just the same as B or delCdt (total uptake by land), this
% variable is used in the diagnostic (single deconv) version of this model
residualLandUptake = [];
residualLandUptake(:,1) = year_ocean2;
dtdelpCO2a = []; %
dtdelpCO2a(:,1) = year_ocean2;


fas = zeros(2516,2);
fas(:,1) = year_ocean2;


for i = 1:length(year_ocean2)-1; % changed this to -1 -- any adverse effects?
    
    % ocean uptake - code from joos_general_fast_annotate2
    
    %Calculate flux 
    % figure out what the fas units are
    fas(i,1) = year_ocean(i);
    % bug: note here fas is just being filled with zeros bc dpco2a and
    % dpco2s are both just full of zeros
    fas(i,2) = (kg/Aoc)*(dpCO2a(i,2) - dpCO2s(i,2)); % air-sea flux of CO2
    %TODO: need to transpose dimensions of dpco2a - not sure where went
    %wrong, check LR code
    % dpco2a comes out of MLOinterp as 2522x2 double
    % for some reason, initializing dpco2a = [1:2522,1] gives me the
    % opposite dimensions I'd expect (I wanted 2521x1, got 1x2521)
    % TODO: multiplying by dpCO2a, but it's just zero at the moment

   % convolve the air-sea flux and the pulse response function (Joos 1996)
    w = conv(fas(1:i,2),r(1:i,2)); 

    % Calculate delDIC 
    % units on this?
    delDIC(i+1,1) = year_ocean(i+1); % filling time column for delDIC
    delDIC(i+1,2) = (c/h)*w(i)*dt; % change in DIC
    %end

    %Calculate dpCO2s from DIC - from Joos 1996
    dpCO2s(i+1,2) = (1.5568 - (1.3993E-2)*T_const)*delDIC(i+1,2) + ...
        (7.4706-0.20207*T_const)*10^(-3)*(delDIC(i+1,2))^2 - ...
        (1.2748-0.12015*T_const)*10^(-5)*(delDIC(i+1,2))^3 + ...
        (2.4491-0.12639*T_const)*10^(-7)*(delDIC(i+1,2))^4 - ...
        (1.5468-0.15326*T_const)*10^(-10)*(delDIC(i+1,2))^5;
    
    
    % land uptake - code from bioboxtwo_sub10_annotate.m for CO2
    
    % fertilization, from bioboxtwo_subN_annotate.m for nitrogen
    % fertilization
    
    % Question/TODO: what are we doing about temp-dependent photosynthesis
    % and respiration? TODO part is to make if cases for these
    % currently nitrogen == 0, so ignore following code Jan 7, 2018
    if nitrogen == 1 % nitrogen fertilization case
        
        % fast box
        % temperature-dependent respiration
        C1dt(i,2) = Ka1*Catm*(1 + eps*dpCO2a(i,2)/Catm + gamma*ff1(i+a-1,2))...
            - K1a*Q1a^((temp_anom(i,2)-T0)/10)*(C1 + delC1(i,2)); 
        % temperature-dependent photosynthesis
        % C1dt(m,2) = Ka1*Catm*(1 + eps*dpCO2a(m,2)/Catm + gamma*ff1(m+a-1,2))...
        %*(1 + Q1a*(temp_anom(m,2)-T0)) - K1a*(C1 + delC1(m,2)); 
        

        % slow box
        % temperature-dependent respiration 
        C2dt(i,2) = Ka2*Catm*(1 + eps*dpCO2a(i,2)/Catm + gamma*ff1(i+a-1,2))...
            - K2a*Q2a^((temp_anom(i,2)-T0)/10)*(C2 + delC2(i,2));
        % temperature-dependent photosynthesis
        % C2dt(m,2) = Ka2*Catm*(1 + eps*dpCO2a(m,2)/Catm + gamma*ff1(m+a-1,2))...
        %*(1 + Q2a*(temp_anom(m,2)-T0)) - K2a*(C2 + delC2(m,2)); 
        
        % box total change in concentrations
        delC1(i+1,2) = sum(C1dt(:,2))*dt;
        delC2(i+1,2) = sum(C2dt(:,2))*dt;

        % total flux into land
        delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
    
    else % co2 fertilization case (nitrogen == 0)

        % land uptake (biobox10)
        % fast box
        % dpco2a = 
        % note: line below throws error "index exceeds matrix dims" in OG
        % driver
        % this leaves C1dt with all zeros - why?
        % temperature-dependent respiration 
        C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2)) - K1a*Q1a^((temp_anom(i,2)-T0)/10)...
            *(C1 + delC1(i,2)); 
        % temperature-dependent photosynthesis
        % C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2))*(1 + Q1a*(temp_anom(i,2)-T0))...
        %- K1a*(C1 + delC1(i,2)); 

        % slow box
        % temperature-dependent respiration 
        C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2)) - K2a*Q2a^((temp_anom(i,2)-T0)/10)...
            *(C2 + delC2(i,2)); 
        % temperature dependent photosynthesis  
        % C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2))*(1 + Q2a*(temp_anom(i,2)-T0))...
        % - K2a*(C2 + delC2(i,2)); 
        
        % box total change in concentrations
        delC1(i+1,2) = sum(C1dt(:,2))*dt;
        delC2(i+1,2) = sum(C2dt(:,2))*dt;

        % total flux into land
        delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
        
    end
    
    B = delCdt; 


if predict == 1 % prognostic case, running forward model
    
%     fas(i,2) = 0;
%     dtdelpCO2a(i,2) = 0;
%     ff1(i,2) = 0;
%     landusemo(i,2) = 0;
    %B(i,2) = 0;
    
     
    
    % dtdelpCO2a: time derivative, not overall increase since preindustrial
    dtdelpCO2a(i,2) =  ff1(i,2) + landusemo(i,2) - Aoc*fas(i,2) - B(i,2); 
    residualLandUptake(i,2) = dtdelpCO2a(i,2) + Aoc*fas(i,2) - ff1(i,2) - landusemo(i,2);
    % new overall increase since preindustrial = previous value + change (dt)
    dpCO2a(i+1,2) = dpCO2a(i,2) + dtdelpCO2a(i,2)/12; 
    


else % predict == 0
    
    % budget, in ppm/yr
    residualLandUptake(i,2) = dtdelpCO2a_obs(i,2) + Aoc*fas(i,2) - ff1(i,2);% - landusemo(i,2); 
    dpCO2a(i+1,2) = dpCO2a_obs(i+1,2); % in ppm
    

end
    
    % with units, if everyting in ppm, make sure fluxes in ppm/yr
    % find where LR integrates dpco2a etc to get absolute atmospheric co2
    % levels - in her vector CO2a, comes from MLOinterp, just from
    % mlospo_meure vector. Only fiddles with dpco2a etc because ocean and
    % land uptake is a function of atmospheric change, not total
    % concentration
    % I can just integrate my change from preinductrial atmospheric co2
    % concentration plus all the changes in dpco2
end 

% some wrap-up from LR's joos_general_fast_annotate2.m which makes it
% easier to compare my variables to the LR variables
zero1 = find(delDIC(:,1) == 0);
delDIC(zero1,2) = NaN;
delDIC(zero1,1) = NaN;

zero2 = find(dpCO2s(:,2) == 0);
dpCO2s(zero2,2)  = NaN;

% calculate the flux for the last time point
fas(length(year_ocean),1) = year_ocean(length(year_ocean));
fas(length(year_ocean),2) = (kg/Aoc)*(dpCO2a(length(year_ocean),2)...
   - dpCO2s(length(year_ocean),2));

% update sign convention
if JLD == 1
    fas(:,2) = -1*fas(:,2);
end


% 10-year smoothing on residualLandUptake
% don't use 10 year mean before 1957, because data are already smoothed
% are start time and end time on l_boxcar indices or years? seems like
% indices, but LR start time of 1225 equates to 1952, not 1957
% somehow this works for LR
i = find(residualLandUptake(:,1) == 1957);
%[residual10a] = l_boxcar(residualLandUptake,10,12,i,length(residualLandUptake),1,2);

%% post-loop processsing and plotting

ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);


if predict == 1
    
    % calculate residualLandUptake with updated value for fas at last point
    residualLandUptake(length(year_ocean),2) = dtdelpCO2a(length(year_ocean),2)...
        + Aoc*fas(length(year_ocean),2) - ff1(length(year_ocean),2);% - landusemo(i,2);

    % mass balance check
    sumCheck = [];
    sumCheck(:,1) = year_ocean2; 
    
    if JLD == 1 % fas is positive into atmosphere
        sumCheck(:,2) = ff1(:,2)  + Aoc*fas(:,2) + residualLandUptake(:,2)...
            - dtdelpCO2a(:,2)+ landusemo(:,2);
    else % fas is positive into ocean
        sumCheck(:,2) = dtdelpCO2a(:,2) - residualLandUptake(:,2) - ff1(:,2)...
            + Aoc*fas(:,2) - landusemo(:,2);
    end
    
    figure('name','sumCheck - PREDICT');
    plot(sumCheck(:,1),sumCheck(:,2));
    axis([1800 2010 -10 10])
    title('sumCheck = atmos - land - ff + Aoc*fas')
    xlabel('Year')
    
    H = figure('name','residualLandUptake plus components JLD - PREDICT');
    if JLD == 1
    plot(residualLandUptake(:,1),residualLandUptake(:,2),'-g',ff1(:,1), ff1(:,2), ...
        '-k', dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r', fas(:,1),Aoc*fas(:,2),'-b');
%     plot(B(:,1),B(:,2),'-g',ff1(:,1), ff1(:,2), ...
%         '-k', dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r', fas(:,1),Aoc*fas(:,2),'-b');
    else
        plot(residualLandUptake(:,1),residualLandUptake(:,2),'-g',ff1(:,1), ff1(:,2), ...
        '-k', dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r', fas(:,1),-Aoc*fas(:,2),'-b');
    end
    axis([1800 2010 -10 10])
    legend('residualLandUptake','fossil fuel','atmosphere','ocean','Location','SouthWest')
    title('residualLandUptake plus components JLD - PREDICT')
    xlabel('Year ')
    ylabel('ppm/year  Positive = source, negative = sink ')


else % predict == 0

    % calculate residualLandUptake with updated value for fas at last point
    residualLandUptake(length(year_ocean),2) = dtdelpCO2a_obs(length(year_ocean),2)...
        + Aoc*fas(length(year_ocean),2) - ff1(length(year_ocean),2);% - landusemo(i,2);

    % mass balance check
    sumCheck = [];
    sumCheck(:,1) = year_ocean2; 
    % the logical version:
    %sumCheck(:,2) = ff1(:,2) + landusemo(:,2) - residualLandUptake(:,2) - Aoc*fas(:,2) - dtdelpCO2a_obs(:,2);
    if JLD == 1 % fas is positive into atmosphere
        sumCheck(:,2) = ff1(:,2)  + Aoc*fas(:,2) + residualLandUptake(:,2) - dtdelpCO2a_obs(:,2); %+ landusemo(:,2)
    else % fas is positive into ocean
        sumCheck(:,2) = dtdelpCO2a_obs(:,2) - residualLandUptake(:,2) - ff1(:,2) + Aoc*fas(:,2);% - landusemo(:,2);
    end
    figure('name','sumCheck - NO PREDICT');
    plot(sumCheck(:,1),sumCheck(:,2));
    axis([1800 2010 -10 10])
    title('sumCheck = atmos - land - ff + Aoc*fas')
    xlabel('Year')
    % yields deviations from 0 on order of 10^-16
    
    
    % plot
    H = figure('name','residualLandUptake plus components JLD - NO PREDICT');
    plot(residualLandUptake(:,1),residualLandUptake(:,2),'-g',ff1(:,1), ff1(:,2), ...
        '-k', dtdelpCO2a_obs(:,1),dtdelpCO2a_obs(:,2),'-r', fas(:,1),-Aoc*fas(:,2),'-b');
    axis([1800 2010 -10 10])
    legend('residualLandUptake','fossil fuel','atmosphere','ocean','Location','SouthWest')
    title('residualLandUptake plus components JLD - NO PREDICT')
    xlabel('Year')
    ylabel('ppm/year  Positive = source, negative = sink')

end

    % arranging location of figures

    vert = 300; %300 vertical pixels
    horz = 600; %600 horizontal pixels
    set(H,'Position',[(3*width/4)-horz/2, (height/2)-vert/2, horz, vert]);
    I = openfig('LR_plot.fig');
    set(I,'Position',[(width/4)-horz/2, (height/2)-vert/2, horz, vert]);
