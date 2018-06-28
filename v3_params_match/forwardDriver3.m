% forward model version 3
%
% creating to match LR fitted params for 2009+7/12
%
% march 21, 2018
% author: julia dohner, with code adapted from lauren rafelski

clear all

%% define time frame, cases

landusedata = 'hough03'; %hough, hansis, hough03, const, const2;
    

varSST = 0; %1 if variable sst, 0 if fixed sst
nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?
filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered
tropicalLU = 1; % 1 = use tropical LU, 0 = extratropical LU
inputStr = 'CHM-V';

Tconst = 18.2; % surface temperature, deg C, from Joos 1996
ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;%2009+(7/12);%
end_year_plot = 2015.5;
year2 = (start_year:(1/ts):end_year)';
beta = [0.5;2]; % initial guesses for model fit (epsilon, q10)
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996


save('yearInfo','start_year','end_year','ts','year2');


%% fitting parameters

%prompt = 'Which case would you like to retrieve fitted parameters for? Please use the form XXX-X with single quotes.  ';
%inputStr = input(prompt);


% high land use: LUlevel = 1
% low land use: LUlevel = 2
%
% high ocean uptake: oceanUptake = 1
% medium ocean uptake: oceanUptake = 2
% low ocean uptake: oceanUptake = 3
%
% temperature dependent: tempDep = 1
% temperature independent: tempDep = 2
if inputStr == 'CHM-C'
    LUlevel = 1;
    oceanUptake = 2;
    tempDep = 0;
    elseif inputStr == 'CHH-V'
        LUlevel = 1;
        oceanUptake = 1;
        tempDep = 1;
    elseif inputStr == 'CLH-V'
        LUlevel = 2;
        oceanUptake = 1;
        tempDep = 1;
    elseif inputStr == 'CHM-V'
        LUlevel = 1;
        oceanUptake = 2; 
        tempDep = 1;
    elseif inputStr == 'CLM-V'
        LUlevel = 2;
        oceanUptake = 2;
        tempDep = 1;
    elseif inputStr == 'CHL-V'
        LUlevel = 1;
        oceanUptake = 3;  
        tempDep = 1;
    elseif inputStr == 'CLL-V'
        LUlevel = 2;
        oceanUptake = 3; 
        tempDep = 1;
    else
        prompt2 = 'Please try again. Which cases would you like to see plotted? Please use the form XXX-X with single quotes.  ';
        inputStr = input(prompt2);
end

if tempDep == 0
    beta = [0.5,1];
end


%% load data

% give access to data files in co2_forward_data folder
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/necessary_data'))';


[dtdelpCO2a_obs,dpCO2a_obs,~,~,CO2a_obs] = getObservedCO2_2(ts,start_year,end_year);


%% get temp record
[temp_anom, ~] = tempRecord2(start_year,end_year,dt);
  
X = temp_anom(:,:);

%% fitting parameters for cases

%[~, ff, LU, LUex] = getSourceSink3(year2, ts); % for LR record
%[ff, LU] = getSourceSink4(year2, ts); % for trying different LU
[ff, LU] = getSourceSink5(year2, ts,landusedata); % for updated FF & LU

if tropicalLU == 0
    LU = LUex;
end

[fas,sstAnom] = jooshildascale_annotate2(start_year,end_year,ts,ff,varSST,Tconst);
%[fas] = joosHILDAvarT_scale_JLD(start_year,end_year,ts,ff,varSST);

% scaling ocean uptake
if oceanUptake == 1 % high ocean uptake
    fas(:,2) = fas(:,2)*1.3;
    disp('ocean multiplied by 1.3')
elseif oceanUptake == 3 % low ocean uptake
        fas(:,2) = fas(:,2)*0.7;
        disp('ocean multiplied by 0.7')
else
    disp('ocean not multiplied')
end

%%----------------------------------------------------------------------%%
% Calculate residual land uptake
% run to 8/2009
% using high land use emissions
% i = length(year2);

i1 = find(dtdelpCO2a_obs(:,1) >= start_year,1);
j1 = find(dtdelpCO2a_obs(:,1) >= end_year);
% i1 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(start_year+(1/24))));
% j1 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(end_year+(1/24))));

% TODO: need to extend dtdelpCO2a record (obs CO2 record) to 2016
%dtdelpCO2a = dtdelpCO2a(i1:j1,:);

i2 = find(ff(:,1) == start_year);
j2 = find(ff(:,1) == end_year);
% i2 = find(dtdelpCO2a(:,1) >= start_year,1);
% j2 = find(dtdelpCO2a(:,1) >= end_year);



%if tropicalLU == 1
% % Extend land use record by making recent emissions equal to last
% % record
% LU(1874:1916,1) = year(1874:1916);
% LU(1874:1916,2) = LU(1873,2);
% % % 

residual(:,1) = year2;
%residual(:,2) = dtdelpCO2a(i1:j1,2) - ff(i2:j2,2)....
%+ Aoc*fas(i3:j3,2) - LU(1:length(year2),2);
residual(:,2) = dtdelpCO2a_obs(:,2) - ff(:,2) + Aoc*fas(:,2) - LU(:,2);

%else 
    
%Extend extratropical emissions by assuming emissions are zero
LUex(1802:length(year2),1) = LU(1802:length(year2),1);
LUex(1802:length(year2),2) = 0;

% using extratropical emissions only
residual2(:,1) = year2;
residual2(:,2) = dtdelpCO2a_obs(:,2) - ff(:,2)....
+ Aoc*fas(:,2) - LUex(:,2);

%end


%%----------------------------------------------------------------------%%

% calculate a 10-year running boxcar average of the residual land uptake
% don't use 10 year mean before 1957, because data are already smoothed
% (ice core)
if(LUlevel==1) %high land use
    %[residual10] = l_boxcar(residual,10,12,1,length(residual),1,2);
    % l_boxcar(func,boxlength,dt,starttime,endtime,datecol,numcol)
    [residual10a] = l_boxcar(residual,10,12,1225,length(residual),1,2);
    residual10(1:1284,:) = residual(1:1284,:);
    residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

    decon = residual;
    
elseif(LUlevel ==2) % low land use
    %[residual10] = l_boxcar(residual2,10,12,1,length(residual2),1,2);
    [residual10a] = l_boxcar(residual2,10,12,1225,length(residual2),1,2);
    residual10(1:1284,:) = residual2(1:1284,:);
    residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

    decon = residual2;
end


%JLD note: residual10(601:end,2) - check dims
%% find model fit using a nonlinear regression

if tempDep == 1
    if(filter == 1) % fit to 10-year filtered record

        [betahat,resid,J] = nlinfit(temp_anom,residual10(601:end,2),'Table2_tempDep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513; change 601 to 1297

    elseif(filter == 2) % fit to unfiltered record

        [betahat,resid,J] = nlinfit(temp_anom,decon(1057:1513,2),'Table2_tempDep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513. After 1958: 1297; was on 1309:end

    end
else
        if(filter == 1) % fit to 10-year filtered record

        [betahat,resid,J] = nlinfit(temp_anom,residual10(601:end,2),'Table2_tempIndep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513; change 601 to 1297

        elseif(filter == 2) % fit to unfiltered record

        [betahat,resid,J] = nlinfit(temp_anom,decon(1057:1513,2),'Table2_tempIndep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513. After 1958: 1297; was on 1309:end

    end
    
end

%% Look at covariances and correlations between model result and calculated land uptake 
[N,P] = size(temp_anom);
% 
covariance = inv(J(1:1177,:)'*J(1:1177,:))*sum(resid(1:1177,:).^2)/(N-P) 

[isize,jsize] = size(covariance);

for k=1:isize
    for j=1:jsize
        correlation(k,j) = covariance(k,j)/sqrt(covariance(k,k)*covariance(j,j));
    end
end

%% Get uncertainties of best fit values
ci = nlparci(betahat,resid,J);

%% Redefine values of epsilon, gamma and Q1


if(nitrogen == 1) % N-fertilization case
    epsilon = 0% betahat(2)
    gamma = betahat(1)
    Q1 = betahat(2)
    Q2 = 1;%betahat(3)
else % co2-fertilization case
    epsilon = betahat(1)%0.79;%
    Q1 = betahat(2)%4.91;
    Q2 = 1;%betahat(2)
end

% epsilon = 0.7743;
% Q1 = 5.8901;


%year2 = year2';

%% plotting // fixing params here
%prompt = 'Plot? Write "yes" or "no" with single quotes.  ';
%inputStr2 = input(prompt);
inputStr2 = 'yes';



% get all records again for new end_year (for full plot)
year3 = (start_year:(1/ts):end_year_plot)';
if end_year ~= end_year_plot
    [dtdelpCO2a_obs,dpCO2a_obs,~,~,CO2a_obs] = getObservedCO2_2(ts,start_year,end_year_plot);
    [temp_anom, ~] = tempRecord2(start_year,end_year_plot,dt);
    [ff, LU] = getSourceSink5(year3, ts,landusedata); % for updated FF & LU
    [fas,sstAnom] = jooshildascale_annotate2(start_year,end_year_plot,ts,ff,varSST,Tconst);
end

% Run the best fit values in the model again to plot


[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year3,dpCO2a_obs,temp_anom); 
 
%%% Nitrogen%%%

% [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN_annotate(epsilon,Q1,Q2,gamma,ff,ts,year2,dpCO2a,X);
   
    
delCdt(:,2) = -delCdt(:,2);

% 10 year moving boxcar average of model result
[delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);

if filter == 1 
    yhat2 = delC10(:,2);
else
    yhat2 = delCdt(:,2);   
end



% Do "reverse deconvolution" to calculate modeled CO2 change
i4 = find(fas(:,1) == start_year);
j4 = find(fas(:,1) == end_year);


if(LUlevel==1)
    newat(:,1) = year3;
    newat(:,2) =  ff(:,2) - Aoc*fas(:,2) + LU(:,2) + delCdt(:,2);

elseif(LUlevel==2)
    
    % 1850 to 2005  --- units on this?
    newat(:,1) = year3;
    newat(:,2) =  ff(:,2) - Aoc*fas(:,2) + LUex(:,2) + delCdt(:,2);

end

    
co2_preind = 600/2.12; % around 283 ppm (preindustrial)

atmcalc(:,1) = year3;
atmcalc(:,2) = cumsum(newat(:,2)/12);
atmcalc(:,2) = atmcalc(:,2)+co2_preind;

co2_diff(:,1) = year3;
co2_diff(:,2) = CO2a_obs(:,2)-atmcalc(:,2);
i6 = find(co2_diff(:,1) == 1959);
j6 = find(co2_diff(:,1) == 1979);
meandiff = mean(co2_diff(i6:j6,2)); % mean difference over 1959-1979
atmcalc2 = atmcalc(:,2)+meandiff;

obsCalcDiff(:,1) = year3;
obsCalcDiff(:,2) = CO2a_obs(:,2) - atmcalc2(:,1); 

if varSST == 1
    
    figure('Name','Modeled vs. Observed CO2')
    subplot(4,1,1)
    plot(CO2a_obs(:,1), CO2a_obs(:,2),year2,atmcalc2);
    xlabel('year')
    ylabel('ppm CO2')
    title('Atmospheric CO2 history')
    legend('Observations','Temperature-dependent model','location','northwest');
    grid
    subplot(4,1,2)
    plot(obsCalcDiff(:,1),obsCalcDiff(:,2));
    line([year3(1),year3(end)],[0,0],'linestyle','--');
    grid
    legend('Observed - modeled CO2','location','northeast')
    xlabel('year')
    ylabel('ppm CO2')
    title('Observed CO2 Deviation from Modeled CO2')
    subplot(4,1,3)
    plot(temp_anom(:,1),temp_anom(:,2))
    line([year3(1),year3(end)],[0,0],'linestyle','--');
    grid
    legend('Temperature anomaly','location','northeast')
    xlabel('year')
    ylabel('deg C')
    title('Temperature anomaly')
    subplot(4,1,4)
    plot(sstAnom(:,1),sstAnom(:,2))
    line([year3(1),year3(end)],[0,0],'linestyle','--');
    grid
    legend('SST anomaly','location','northeast')
    xlabel('year')
    ylabel('deg C')
    title('SST anomaly')
    
else
    

    figure('Name','Modeled vs. Observed CO2')
    subplot(3,1,1)
    plot(CO2a_obs(:,1), CO2a_obs(:,2),year3,atmcalc2);
    xlabel('year')
    ylabel('ppm CO2')
    title('Atmospheric CO2 history')
    legend('Observations','Temperature-dependent model','location','northwest');
    grid
    subplot(3,1,2)
    plot(obsCalcDiff(:,1),obsCalcDiff(:,2));
    line([year3(1),year3(end)],[0,0],'linestyle','--');
    grid
    legend('Observed - modeled CO2','location','northeast')
    xlabel('year')
    ylabel('ppm CO2')
    title('Observed CO2 Deviation from Modeled CO2')
    subplot(3,1,3)
    plot(temp_anom(:,1),temp_anom(:,2))
    grid
    legend('Temperature anomaly','location','northeast')
    xlabel('year')
    ylabel('deg C')
    title('Temperature anomaly')

end
saveas(gcf,'CO2recordsFig.fig') 

figure('Name','Land Uptake')
plot(residual10(:,1),residual10(:,2),delC10(:,1),yhat2,LU(:,1), LU(:,2),'-.')
line([year3(1),year3(end)],[0,0],'linestyle',':');
set(gca,'Xlim',[1850 2010]) 
title('Residual Land Flux')
legend('Residual land flux','Residual land flux (model)','Land use emissions','location','northwest')
xlabel('year')
ylabel('ppm / year')
grid

saveas(gcf,'landFluxFig.fig')

atmos_smooth = smooth(CO2a_obs(:,2),5);
dtdelpCO2a_sm = smooth(dtdelpCO2a_obs(:,2),29); %~ 5 yr running mean. boxcar of 59 months
imbalance = 0;

% sources and sinks plot in fluxes
figure('Name','Sources and Sinks')
plot(ff(:,1),ff(:,2),year3,-dtdelpCO2a_sm,fas(:,1),...
    -Aoc*fas(:,2),delC10(:,1),yhat2,LU(:,1),LU(:,2))
line([year3(1),year3(end)],[0,0],'linestyle',':');
set(gca,'Xlim',[1900 2010]) 
legend('fossil fuel','observed atmosphere','ocean','modeled land','land use','Location','SouthWest')

if end_year == 2015.5 && end_year_plot == 2015.5 && strcmp(landusedata,'const2')
    if tempDep == 1 
            Vpresent_constant2_co2 = atmcalc2;
            Vpresent_constant2_resid = obsCalcDiff;
            Vpresent_constant2_q10 = Q1;
            Vpresent_constant2_eps = epsilon;
            Vpresent_constant2_year = year2;
            save('timeframe_constant2_presentV','Vpresent_constant2_co2','Vpresent_constant2_resid',...
                'Vpresent_constant2_q10','Vpresent_constant2_eps','Vpresent_constant2_year')
    elseif tempDep == 0
            Cpresent_constant2_co2 = atmcalc2;
            Cpresent_constant2_resid = obsCalcDiff;
            Cpresent_constant2_q10 = Q1;
            Cpresent_constant2_eps = epsilon;
            Cpresent_constant2_year = year2;
            save('timeframe_constant2_presentC','Cpresent_constant2_co2','Cpresent_constant2_resid',...
                'Cpresent_constant2_q10','Cpresent_constant2_eps','Cpresent_constant2_year')
    end

end
        %%
        


if end_year == 2005.5 && end_year_plot == 2015.5
    if tempDep == 1
        
        if strcmp(landusedata,'hough')
            V2005_hough_co2 = atmcalc2;
            V2005_hough_resid = obsCalcDiff;
            V2005_hough_q10 = Q1;
            V2005_hough_eps = epsilon;
            V2005_hough_year = year2;
            save('timeframe_hough_2005V','V2005_hough_co2','V2005_hough_resid',...
                'V2005_hough_q10','V2005_hough_eps','V2005_hough_year')
            
        elseif strcmp(landusedata,'hansis')
            V2005_hansis_co2 = atmcalc2;
            V2005_hansis_resid = obsCalcDiff;
            V2005_hansis_q10 = Q1;
            V2005_hansis_eps = epsilon;
            V2005_hansis_year = year2;
            save('timeframe_hansis_2005V','V2005_hansis_co2','V2005_hansis_resid',...
                'V2005_hansis_q10','V2005_hansis_eps','V2005_hansis_year')

        elseif strcmp(landusedata,'hough03')
            V2005_hough03_co2 = atmcalc2;
            V2005_hough03_resid = obsCalcDiff;
            V2005_hough03_q10 = Q1;
            V2005_hough03_eps = epsilon;
            V2005_hough03_year = year2;
            save('timeframe_hough03_2005V','V2005_hough03_co2','V2005_hough03_resid',...
                'V2005_hough03_q10','V2005_hough03_eps','V2005_hough03_year')
      
        elseif strcmp(landusedata,'const')
            V2005_constant_co2 = atmcalc2;
            V2005_constant_resid = obsCalcDiff;
            V2005_constant_q10 = Q1;
            V2005_constant_eps = epsilon;
            V2005_constant_year = year2;
            save('timeframe_constant_2005V','V2005_constant_co2','V2005_constant_resid',...
                'V2005_constant_q10','V2005_constant_eps','V2005_constant_year')
        
        elseif strcmp(landusedata,'const2')
            V2005_constant2_co2 = atmcalc2;
            V2005_constant2_resid = obsCalcDiff;
            V2005_constant2_q10 = Q1;
            V2005_constant2_eps = epsilon;
            V2005_constant2_year = year2;
            save('timeframe_constant2_2005V','V2005_constant2_co2','V2005_constant2_resid',...
                'V2005_constant2_q10','V2005_constant2_eps','V2005_constant2_year')
   
        
        end
    elseif tempDep == 0
               
        if strcmp(landusedata,'hough')
            C2005_hough_co2 = atmcalc2;
            C2005_hough_resid = obsCalcDiff;
            C2005_hough_q10 = Q1;
            C2005_hough_eps = epsilon;
            C2005_hough_year = year2;
            save('timeframe_hough_2005C','C2005_hough_co2','C2005_hough_resid',...
                'C2005_hough_q10','C2005_hough_eps','C2005_hough_year')
            
        elseif strcmp(landusedata,'hansis')
            C2005_hansis_co2 = atmcalc2;
            C2005_hansis_resid = obsCalcDiff;
            C2005_hansis_q10 = Q1;
            C2005_hansis_eps = epsilon;
            C2005_hansis_year = year2;
            save('timeframe_hansis_2005C','C2005_hansis_co2','C2005_hansis_resid',...
                'C2005_hansis_q10','C2005_hansis_eps','C2005_hansis_year')
 
        elseif strcmp(landusedata,'hough03')
            C2005_hough03_co2 = atmcalc2;
            C2005_hough03_resid = obsCalcDiff;
            C2005_hough03_q10 = Q1;
            C2005_hough03_eps = epsilon;
            C2005_hough03_year = year2;
            save('timeframe_hough03_2005C','C2005_hough03_co2','C2005_hough03_resid',...
                'C2005_hough03_q10','C2005_hough03_eps','C2005_hough03_year')
      
        elseif strcmp(landusedata,'const')
            C2005_constant_co2 = atmcalc2;
            C2005_constant_resid = obsCalcDiff;
            C2005_constant_q10 = Q1;
            C2005_constant_eps = epsilon;
            C2005_constant_year = year2;
            save('timeframe_constant_2005C','C2005_constant_co2','C2005_constant_resid',...
                'C2005_constant_q10','C2005_constant_eps','C2005_constant_year')
            
        elseif strcmp(landusedata,'const2')
            C2005_constant2_co2 = atmcalc2;
            C2005_constant2_resid = obsCalcDiff;
            C2005_constant2_q10 = Q1;
            C2005_constant2_eps = epsilon;
            C2005_constant2_year = year2;
            save('timeframe_constant2_2005C','C2005_constant2_co2','C2005_constant2_resid',...
                'C2005_constant2_q10','C2005_constant2_eps','C2005_constant2_year')
        end

    end
end

% 
% 
% if strcmp(landusedata,'hough')
%     if tempDep == 1
%         %if end_year == 2015.5
%         %    save('timeframe_presentV','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
%         if end_year == 2005.5
%             V2005_hough_co2 = atmcalc2;
%             V2005_hough_resid = obsCalcDiff;
%             V2005_hough_q10 = Q1;
%             V2005_hough_eps = epsilon;
%             V2005_hough_year = year2;
%             save('timeframe_hough_2005V','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
% %         elseif end_year == 1995.5
% %             save('timeframe_1995V','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
% %         elseif end_year == 1975.5
% %             save('timeframe_1975V','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
%         end
%     elseif tempDep == 0
%         if end_year == 2015.5
%             save('timeframe_presentC','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
%         elseif end_year == 2005.5
%             save('timeframe_2005C','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
%         elseif end_year == 1995.5
%             save('timeframe_1995C','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
%         elseif end_year == 1975.5
%             save('timeframe_1975C','obsCalcDiff','atmcalc2','Q1','epsilon','year2')
%         end
%     end
%         
% end

    
