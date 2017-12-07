% file ForwardModel_Driver2.m
% version where data starts at 1850, ends at 2006 across all files

clear all; %close all;

% give access to data files in co2_forward_data folder
addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));

% predict = 1 --> prognostic, calculating dpCO2a in motherloop
% predict = 0 --> diagnostic (single deconvolution), feeding dpCO2a from
% MLOinterp, calculating residual land uptake in motherloop
predict = 0;


%% set up ocean

% ~~~~~~~copying from jooshildascale_annotate2.m from here down~~~~~~~

ts = 12; % number of data points/year, 12 in both land and ocean
start_year_ocean = 1800;
end_year_ocean = 2006; % this is the end of landusemo (tland ends at 2006+10/12, extratrop ends at 2000)
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas ex?change rate, yr^-1, from Joos 1996

% don't get atmospheric co2 record from MLOinterp
[ff1] = load_fossil2(ts); % get fossil fuel emissions

% initializing year vector (from MLOinterpolate_increment2.m)
dt = 1/ts;
year_ocean = start_year_ocean:dt:end_year_ocean;
year_ocean = year_ocean - 0.375;
% adjusting year vector to match year vector from LR ocean code

% back to jooshildascale_annotate2.m

% Response function to calculate ocean uptake
[t,r] = HILDAResponse(year_ocean);

% in the original code, would call joos_general_fast_annotate2 here, but
% I'm putting part of that function into the motherloop, so here I'm just
% calling the prep lines from that function (pre-loop)

%%%%%%% next step:

% dpCO2a is initialized in MLOinterp, but we need it for its dimensions to
% initialize dpCO2s
% dpCO2a as outputted from MLOinterp is a 2522x2 double
% (just startyear:1/12:end year)
%dpCO2a = zeros(2521,2); %[1:2522,2];

% this is how it's initialized in MLOinterp:
% % Calculate change in atmospheric concentration
% find floor of 100*all rows in first column of MLOice (so time?) - want to
% find where to start using the record, which is 1850.04 (floor(100*start
% year)
% MLOiceinterp is all values from CO22 = interp1(years0,CO2,years,'spline');
% interpolating from years values provided with mlosp_meure (OG data) to
% new years vector, so CO22 is just this interpolated data

% CO2 data start = 1640
% end = 2010.1
% only use data from 1850 onwards, thus 1922 long
% but we have data from before 1850- in fact, all the way back to 1640
% maybe use 1640 to 1800 data as starter values for dpco2a using same
% interpolation method as MLO interp

% dpCO2a is the change in atmospheric CO2 from preindustrial value
dpCO2a = zeros(length(year_ocean),2); %[1:2483,2];
dpCO2a(:,1) = year_ocean; % fill the time column

dpCO2a_LRdata = load('dpCO2a_LRocean.mat');
dpCO2s_LRdata = load('dpCO2s_LRocean.mat');
dpCO2a_LRocean = dpCO2a_LRdata.dpCO2a;
dpCO2s_LRocean = dpCO2s_LRdata.dpCO2s;

dpCO2a(1:10,2) = dpCO2a_LRocean(1:10,2);

% Dec 5 2017: deleted all the interpolation stuff when trying to figure out a starter
% value


dpCO2s = zeros(length(dpCO2a),2); % dissolved CO2
dpCO2s(:,1) = dpCO2a(:,1);
% dpCO2s(1,2) = 280; % setting initial value (in ppm)
%fas = zeros(length(year),2); % removed because is output by get source
%sink
integral = zeros(length(year_ocean),length(year_ocean));
delDIC = zeros(length(year_ocean),2);

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

[landusemo,ff1,fas,extratrop_landmo] = getsourcesink_scale4; % changed from gss3 to gss4 

%ts = 12; % took out because already defined in ocean
%start_year_land = 1850; % change to start_year_land
%end_year_land = 2009+(7/12); % change to end_year_land
start_year_land = start_year_ocean; % change to start_year_land
end_year_land = end_year_ocean; % change to end_year_land


beta = [0.5;2]; % initial guesses for model fit
% leaving out the betahat stuff because that's just output from the model
% fit using the nonlin regression, but we're already using the best fit
% values here

% leaving out a call to MLO interpolate to get dtdelpCO2a,dpCO2a,year,dt,CO2a
% we get a year vector from MLOinterp, but same as line 26 so won't repeat
% oct 28 2017: changing to add a land year vector for the meantime since
% the land and ocean year vectors are diffferent, and land year vector is
% called in landusemo calls
year_land = start_year_land:dt:end_year_land;

% TODO: the following (extending records) can all happen in getsourcesink_scale3, but wait to
% change until can get the code running and working

% Extend land use record by making recent emissions equal to last
% record
% Question: where do these indices come from???
% if this is counting from 1850, then the index 1874 corresponds to year 2006.2
% if this is counting from 1850, then the index 1916 corresponds to year 2009.7
% if this is counting from 1850, then the index 1873 corresponds to year 2006.1

% if this is counting from 1800, then the index 1916 corresponds to year 1959.7
% if this is counting from 1800, then the index 1874 corresponds to year 1956.2
% if this is counting from 1800, then the index 1916 corresponds to year 1956.1
% most probable: counting from 1850, so tacking on values for 2006.2 onto 2009.7
% so this is dealing with stuff after 2006, and seems unnecessary
% landusemo(1874:1916,1) = year_land(1874:1916); % TODO: this is calling on land_year vector, but right now calls ocean_year
% landusemo(1874:1916,2) = landusemo(1873,2);

% Question: why didn't LR just extend the interpolation for landusemo and
% extratop? Rather than what she did where she set the recent values equal
% to something

% % 
%Extend extratropical emissions by assuming emissions are zero
% extratrop_landmo(1802:1916,1) = landusemo(1802:1916,1);
% extratrop_landmo(1802:1916,2) = 0;

% left out: Calculate residual land uptake section, fine to leave out

%%%%%%%%% Nov 8 stopping here %%%%%%%%%

% tland4: the temperature record started at 1880. tland4 extends the
% record back to 1800 by using the mean temperature for the first year
% 
% do a moving boxcar average of the land temperature: 1 year average
% note: in this case the box length (1; second term in l_boxcar) is in
% units of years. dt (12, third term) is the number of points per year
% first column of avg_temp gives the date, second column gives the moving
% average of the land temperature
[avg_temp] = l_boxcar(tland4,1,12,1,2483,1,2); 

avg_temp(1:6,2) = avg_temp(7,2); % make the first 6 points 
% Question: why does avg_temp has zeros for first 6 time points?
% do these time points just not get called? only the corresponding values?

% 10 year moving boxcar average of land temperature
% [avg_temp] = l_boxcar(tland4,10,12,1,2483,1,2);
% 
% avg_temp(1:60,2) = avg_temp(61,2);

%----------------------------------------------------------------------%%
% Pick the temperature record to use
%----------------------------------------------------------------------%%
% 
% ***Use these for various T tests***

% JLD code
% 
% % values for year 1800 to 1850.5 in temp_anom should all be the first value
% % of the avg_temp vector
% % 50*12 = 600 + 6 (for the half year) = first 606 values
% tacking 600 points onto all these calculations
temp_anom(1:606,1) = avg_temp(1:606,1); % fill time column of temp_anom
temp_anom(1:606,2) = landtglob(1,2); %landtglob are temp anomalies
% % TODO: what about filling entire time column?
% % temp_anom(:,1) = avg_temp(:,1); % filling entire time column
% % landtglob starts at year 1850.5, goes until 2010 + 5/12
% accounted for landtglob starting at 1850.5 (rather than 1850) by filling
% in temp_anom starting at 607, rather than 601
temp_anom(607:2516,1) = landtglob(1:1910,1);
temp_anom(607:2516,2) = landtglob(1:1910,1);


%%%%%% LR code %%%%%%

% % first column of temp_anom is just time values
%  temp_anom(1:6,1) =  avg_temp(601:606,1); %Jan 1850-May 1850
%  % want to figure out why put first value in landtglob as first 6 values in
%  % temp anom
%  temp_anom(1:6,2) = landtglob(1,2); %355 instead of 1, 360 instead of 6
%  
%  % here she's dealing with the fact that landtglob starts at 1850.5,
%  % filling this back by 6 months
%  % Ralph: get rid of the lines of code dealing with first 6 months, instead
%  % fill it all the way back to 1800 (in landtglob), so don't mess with here
%  % update T0 to be intiial values from 1800 to 1850, but make sure not same
%  % value that she used for T0
% % update: change her code so it fills all the way back to 1850, have
% % 50*12=726 additional values
%  temp_anom(7:1916,1) = landtglob(1:1910,1); % Starts at the year 1850.5. 
%  temp_anom(7:1916,2) = landtglob(1:1910,2); % 
 
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:6,2) = npptglob(1,2); %355 instead of 1, 360 instead of 6
%  
%  temp_anom(7:1867,1) = npptglob(1:1861,1); % Starts at the year 1850.5
%  temp_anom(7:1867,2) = npptglob(1:1861,2); % 361 instead of 7, 355 instead of 1
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:6,2) = npptno(1,2);
% 
%  temp_anom(7:1867,1) = npptno(1:1861,1); % Starts at the year 1850.5
%  temp_anom(7:1867,2) = npptno(1:1861,2);
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:15,2) = npptso(10,2);
% 
%  temp_anom(7:1867,1) = npptso(1:1861,1); % Starts at the year 1850.5
%  temp_anom(16:1867,2) = npptso(10:1861,2);
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:18,2) = nppttro(13,2);
% 
%  temp_anom(7:1867,1) = nppttro(1:1861,1); % Starts at the year 1850.5
%  temp_anom(19:1867,2) = nppttro(13:1861,2);

%  temp_anom(:,1) = T(601:2467,1);
%  temp_anom(:,2) = T(601:2467,2);
%  
%X = temp_anom(:,:); % just going to call on temp_anom in biobox loop
%portion (rather than calling it as T or X (defined as T, but fed in as X
%in LR code)
 

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

%% OBSERVATION DIAGNOSTIC DEBUGGING SECTION

% MLOinterp isn't called anywhere else in the code
[dtdelpCO2a_obs,dpCO2a_obs,year_obs,dt_obs,CO2a_obs] = MLOinterpolate_increment2(ts,start_year_ocean,end_year_ocean); 

% truncate obs output to match my timeframe
%i = find(floor(100*MLOSPOiceinterp(:,1)) == floor(100*(start_year+(1/24))));
% startIndex = find(floor(dtdelpCO2a_obs(:,1) == start_year_ocean);
% endIndex = find(dtdelpCO2a_obs(:,1) == end_year_ocean);
% dtdelpCO2a_obs = dtdelpCO2a_obs(startIndex:endIndex,:);
% startIndex = 1919; % this is the closest date to 1800- it's actually halfway through december of 1799

% this is bootleg, but it'll have to do for now:
dtdelpCO2a_obs = dtdelpCO2a_obs(1919:4391,:); % indices of closest values to 1800 and 2006


% % I want to interpolate this so that it matches my year vector....
% years0 = dtdelpCO2a_obs(:,1);
% years0 = years0(7:length(years0));
% dtdelpCO2a_obsVals = dtdelpCO2a_obs(:,2);
% dtdelpCO2a_obsVals = dtdelpCO2a_obsVals(7:length(dtdelpCO2a_obsVals));
% % want to interpolate dtdelpCO2a_obs to match my year vector
% dtdelpCO2a_interp = interp1(years0,dtdelpCO2a_obsVals,year_ocean,'spline');
% 
% dtdelpCO2a_obsInterp(:,1) = year_ocean;
% dtdelpCO2a_obsInterp(:,2) = dtdelpCO2a_interp(:,1);




%% probably time for THE MOTHERLOOP

%before this call, year is a 1x2521 double vector (one row vector)
year_ocean2 = year_ocean';


% everything below is just the loop contents from
% bioboxtwo_sub10_annotate.m and joos_general_fast_annotate2, but lacks the
% stuff from both of those functions that help set up the loop.

% i loops through 2521 values

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

if nitrogen == 1 % nitrogen fertilization case, set epsilon = 0 to not allow CO2 fertilization
    eps = 0;
    gamma = beta(1); % gamma is an input only for bioboxtwo_subN (N fertilization)
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
delC1(length(year_land_trans)+1,1) = year(length(year_land_trans),1)+dt;
delC2(length(year_land_trans)+1,1) = year(length(year_land_trans),1)+dt;



% this is just the same as B or delCdt (total uptake by land), this
% variable is used in the diagnostic (single deconv) version of this model
residualLandUptake = [];
residualLandUptake(:,1) = year_ocean2;
dtdeldpCO2a = []; %



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

   
    w = conv(fas(1:i,2),r(1:i,2)); % convolve the air-sea flux and the pulse response function, as in Joos 1996

    % Calculate delDIC 
    delDIC(i+1,1) = year_ocean(i+1); % filling time column for delDIC
    delDIC(i+1,2) = (c/h)*w(i)*dt; % change in DIC
    %end

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
        C1dt(i,2) = Ka1*Catm*(1 + eps*dpCO2a(i,2)/Catm + gamma*ff1(i+a-1,2)) - K1a*Q1a^((temp_anom(i,2)-T0)/10)*(C1 + delC1(i,2)); % temperature-dependent respiration
        % C1dt(m,2) = Ka1*Catm*(1 + eps*dpCO2a(m,2)/Catm + gamma*ff1(m+a-1,2))*(1 + Q1a*(temp_anom(m,2)-T0)) - K1a*(C1 + delC1(m,2)); % temperature-dependent photosynthesis

        % slow box
        C2dt(i,2) = Ka2*Catm*(1 + eps*dpCO2a(i,2)/Catm + gamma*ff1(i+a-1,2)) - K2a*Q2a^((temp_anom(i,2)-T0)/10)*(C2 + delC2(i,2)); % temperature-dependent respiration 
        % C2dt(m,2) = Ka2*Catm*(1 + eps*dpCO2a(m,2)/Catm + gamma*ff1(m+a-1,2))*(1 + Q2a*(temp_anom(m,2)-T0)) - K2a*(C2 + delC2(m,2)); % temperature-dependent photosynthesis

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
        C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2)) - K1a*Q1a^((temp_anom(i,2)-T0)/10)*(C1 + delC1(i,2)); % temperature-dependent respiration 
        % C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2))*(1 + Q1a*(temp_anom(i,2)-T0)) - K1a*(C1 + delC1(i,2)); % temperature-dependent photosynthesis

        % slow box
        C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2)) - K2a*Q2a^((temp_anom(i,2)-T0)/10)*(C2 + delC2(i,2)); % temperature-dependent respiration 
        % C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2))*(1 + Q2a*(temp_anom(i,2)-T0)) - K2a*(C2 + delC2(i,2)); % temperature dependent photosynthesis  

        % box total change in concentrations
        delC1(i+1,2) = sum(C1dt(:,2))*dt;
        delC2(i+1,2) = sum(C2dt(:,2))*dt;

        % total flux into land
        delCdt(i,2) = C2dt(i,2) + C1dt(i,2);
        
    end
    
    B = delCdt; 


if predict == 1 % prognostic case, running forward model
    
    dtdeldpCO2a(i,2) =  ff1(i,2) + landusemo(i,2) - Aoc*fas(i,2) - B(i,2); %time derivative, not overall increase since preindustrial
    dpCO2a(i+1,2) = dpCO2a(i,2) + dtdeldpCO2a(i,2)/12; % new overall increase since preindustrial = previous value + change (dt)

else % predict == 0
    
    % deconvolution calculation (same as in LR code), should yield same
    % result
    
    % 
    residualLandUptake(i,2) = dtdelpCO2a_obs(i,2) + Aoc*fas(i,2) - ff1(i,2);%  - landusemo(i,2); % budget, in ppm/yr
    % note: LR's dpco2a starts a little before 1800, so visdiff yields
    % different vectors
    % her dpco2a dates come from the find with mlomeure_spo in MLOinterp 
    dpCO2a(i+1,2) = dpCO2a_obs(i+1,2); % in ppm

end
    
    %%%%%%%dpco2a_dt(i,2) =  ff1 + landuse - aoc*fas - modeled land sink %time derivative, not overall increase since preindustrial
% annual increase in ppm/yr (check), calculated monthly (at monthly
% resolution) - therefore increment at current timestep is 1/12 of this -
% RFK

%%%%%%%dpCO2a(i+1,2) = dpCO2a(i,2) + dpCO2a_dt/12;


    % with units, if everyting in ppm, make sure fluxes in ppm/yr
    % find where LR integrates dpco2a etc to get absolute atmospheric co2
    % levels - in her vector CO2a, comes from MLOinterp, just from
    % mlospo_meure vector. Only fiddles with dpco2a etc because ocean and
    % land uptake is a function of atmospheric change, not total
    % concentration
    % I can just integrate my change from preinductrial atmospheric co2
    % concentration plus all the changes in dpco2
end 

% debug plots
% include LR plots
% code to always pull up graphs from both
figure('name','residualLandUptake plus components JLD')

plot(residualLandUptake(:,1),residualLandUptake(:,2),'g',ff1(:,1),ff1(:,2),'k',dtdelpCO2a_obs(:,1),dtdelpCO2a_obs(:,2),'r',fas(:,1),-Aoc*fas(:,2),'b');
axis([1800 2010 -10 10])
legend('residualLandUptake','fossil fuel','atmosphere','ocean','Location','SouthWest')
title('residualLandUptake plus components JLD ')
xlabel('Year ')
ylabel('ppm/year  Positive = source, negative = sink ')
openfig('LR_plot.fig');
%ff1(:,1),ff1(:,2),'-k',dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r',fas(:,1),-Aoc*fas(:,2),'-b',landflux(:,1),landflux(:,2),'-g',year(1,:),x,'--k')

%residualLandUptake(i,2) = dtdelpCO2a_obs(i,2) - ff1(i,2) + Aoc*fas(1,2) - landusemo(1,2); % budget, in ppm/yr
%plot(ff1(:,1),ff1(:,2),'-k',dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r',fas(:,1),-Aoc*fas(:,2),'-b',landflux(:,1),landflux(:,2),'-g',year(1,:),x,'--k')

% post-loop code from joos_general_fast_annotate2.m
% calculate the flux for the last time point
fas(length(year_ocean),1) = year_ocean(length(year_ocean));
fas(length(year_ocean),2) = (kg/Aoc)*(dpCO2a(length(year_ocean),2) - dpCO2s(length(year_ocean),2));

%plot(dpCO2s(:,1),dpCO2s(:,2));

% zero1 = find(delDIC(:,1) == 0);
% delDIC(zero1,2) = NaN;
% delDIC(zero1,1) = NaN;
% 
% zero2 = find(dpCO2s(:,2) == 0);
% dpCO2s(zero2,2)  = NaN;


%% Misc from nonlin_land_Qs_annotate.m

% Leaving all of this out. Question: probably don't need, right?
% % Run the best fit values (parameters from Table 2) in the model again to plot
% if nitrogen == 1
%     [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN_annotate(epsilon,Q1,Q2,gamma,ff1,ts,year2,dpCO2a,temp_anom);
% else 
%     [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,temp_anom); 
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




