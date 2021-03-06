% forward model version 3
%
% creating to match LR fitted params for 2009+7/12
%
% march 21, 2018
% author: julia dohner, with code adapted from lauren rafelski

clear all

%% fitting parameters

prompt = 'Which case would you like to retrieve fitted parameters for? Please use the form XXX-X with single quotes.  ';
inputStr = input(prompt);
% high land use: LUlevel = 1
% low land use: LUlevel = 2
%
% high ocean uptake: oceanUptake = 1
% medium ocean uptake: oceanUptake = 2
% low ocean uptake: oceanUptake = 3
%
% temperature dependent: tempDep = 1
% temperature independent: tempDep = 2
if inputStr == 'CHH-V'
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

%% define time frame, cases

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2009+(7/12);%2015.5;% % latest can go is 2015.5 (6 mo before 2016)
year2 = (start_year:(1/ts):end_year)';

% ocean constants
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996

% define cases (nitrogen, filter, tropical vs extratrop lu)
nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?
filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered
tropicalLU = 1; % 1 = use tropical LU, 0 = extratropical LU

beta = [0.5;2]; % initial guesses for model fit (epsilon, q10)
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996

save('caseDetails','ts','start_year','end_year','year2','nitrogen')
% % Initial conditions - set according to Joos et al
% 
% ts = 12;
% start_year = 1850;
% end_year = 2015.5;%2009+(7/12);
% year2 = (start_year:(1/ts):end_year)';
% 
% CO2 = 1; % 1 = CO2 fertilization, 0 = N fertilization model


%% load data

% give access to data files in co2_forward_data folder
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data'));
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data_2016'));



% loading temp data
%addpath(genpath('/Users/juliadohner/Documents/MATLAB/Rafelski_LandOceanModel_JLDedits/co2_forward/co2_forward_data'));


%load land_temp.mat % land temperature records
%load npp_T.mat % NPP-weighted temperature record
%load landwt_T_2011.mat % land temperature anomaly


%[dtdelpCO2a,dpCO2a,~,~,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year);
[dtdelpCO2a,dpCO2a,CO2a,co2_preind] = MLOinterpolate_increment2_recent(ts,start_year,end_year); 



%% get temp record

% X is temp_anom record
[temp_anom, ~] = tempRecord2(start_year,end_year,dt);

 %% fitting parameters for cases

[ff, LU, LUex] = getSourceSink3(year2, ts);

[fas,dpCO2s] = joos_oceanuptake(year2,c,h,kg,T,Aoc,dt);
load fas_LR.mat;
fas2 = fas;

% scaling ocean uptake
if oceanUptake == 1 % high ocean uptake
    fas(:,2) = fas2(:,2)*1.3;
    disp('ocean multiplied by 1.3')
elseif oceanUptake == 3 % low ocean uptake
        fas(:,2) = fas2(:,2)*0.7;
        disp('ocean multiplied by 0.7')
else
    disp('ocean not multiplied')
end

%%----------------------------------------------------------------------%%
% Calculate residual land uptake
% run to 8/2009
% using high land use emissions
% i = length(year2);


%i1 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(start_year+(1/24))));
%j1 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(end_year+(1/24))));
% TODO: need to extend dtdelpCO2a record (obs CO2 record) to 2016
%dtdelpCO2a = dtdelpCO2a(i1:j1,:);
i1 = find(dtdelpCO2a(:,1) == start_year);
j1 = find(dtdelpCO2a(:,1) == end_year);

i2 = find(ff(:,1) == start_year);
j2 = find(ff(:,1) == end_year);

i3 = find(fas(:,1) == start_year);
j3 = find(fas(:,1) == end_year);


%if tropicalLU == 1
% % Extend land use record by making recent emissions equal to last
% % record
% LU(1874:1916,1) = year(1874:1916);
% LU(1874:1916,2) = LU(1873,2);
% % % 

residual(:,1) = year2;
residual(:,2) = dtdelpCO2a(i1:j1,2) - ff(i2:j2,2)....
+ Aoc*fas(i3:j3,2) - LU(1:length(year2),2);

%else 
    
%Extend extratropical emissions by assuming emissions are zero
LUex(1802:length(year2),1) = LU(1802:length(year2),1);
LUex(1802:length(year2),2) = 0;

% using extratropical emissions only
residual2(:,1) = year2;
residual2(:,2) = dtdelpCO2a(i1:j1,2) - ff(i2:j2,2)....
+ Aoc*fas(i3:j3,2) - LUex(1:length(year2),2);

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

%%  PLOTTING WITH FITTED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redefine values of epsilon, gamma and Q1


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


%year2 = year2';

%% plotting 
prompt = 'Plot? Write "yes" or "no" with single quotes.  ';
inputStr2 = input(prompt);

if strcmpi('yes',inputStr2)
% Run the best fit values in the model again to plot


[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,temp_anom); 
 
%%% Nitrogen%%%
if nitrogen == 1
    [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN_annotate(epsilon,Q1,Q2,gamma,ff,ts,year2,dpCO2a,temp_anom);
end
    
delCdt(:,2) = -delCdt(:,2);

% 10 year moving boxcar average of model result
[delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);

% Quantify how good the fit is using the mean squared error (MSE)       
% Use filtered or unfiltered yhat
%----------------------------------------------------------------%
%
% Filtered
%
%------------------------------------------------

if(filter == 1)

yhat2 = delC10(:,2);

%yhat2 = yhat2 + (0 - yhat2(1));
% JLD code below:
% % i5 = find(year2 == 1855);
% %     
% %     e = delC10(i5:end,2) - residual10(i5:end,2); % look at MSE for 1855 to 2000
% %     
% %     misfit = e'*e/length(delC10(i5:end,2));  
% %     
% %     % 1900 to 2005
% %     i6 = find(year2 == 1900);
% %     e2 = delC10(i6:end,2) - residual10(i6:end,2); % look at MSE for 1900 to 2000
% %     misfit2 = e2'*e2/length(delC10(i6:end,2));  
% % 
% % 
% %    error1 = betahat(1)-ci(1);

e = delC10(61:1806,2) - residual10(61:1806,2); % look at MSE for 1855 to 2000

misfit = e'*e/length(delC10(61:1806,2));  

e2 = delC10(601:1806,2) - residual10(601:1806,2); % look at MSE for 1900 to 2000
misfit2 = e2'*e2/length(delC10(601:1806,2));  


error1 = betahat(1)-ci(1);
   
%   error2= betahat(2)-ci(2)

%   error3= betahat(3)-ci(3)
%    
%   figure
% plot(decon(:,1),decon(:,2),delCdt(:,1),yhat2)
% xlabel('year')
% ylabel('ppm CO2/year')
% title('land uptake')
% legend('Residual uptake','land uptake without T effects','land uptake with T effects')

   
    C = cov(residual10(601:(end-1),2),yhat2(601:(end-1),1));
%    
   [R,P,RLO,RUP] = corrcoef(yhat2(601:(end-1),1),residual10(601:(end-1),2));
%   
   R(1,2)^2;

figure
plot(residual10(:,1),residual10(:,2),delC10(:,1),yhat2)
xlabel('year')
ylabel('ppm CO2/year')
title('land uptake')
legend('Residual uptake','land uptake with T effects','location','northwest')
set(gca,'Xlim',[1850 2010])  
grid

%----------------------------------------------------------------%
%
% Unfiltered
%
%----------------------------------------------------------------%
elseif(filter == 2)

yhat2 = delCdt(:,2);    
  
      e = delCdt(61:1806,2) - decon(61:1806,2); % look at MSE for 1855-2000
    
   % misfit = e'*e/length(delCdt(61:1806,2))  
    
    e2 = delCdt(601:1806,2) - decon(601:1806,2); % look at MSE for 1900-2000
    
   % misfit2 = e2'*e2/length(delCdt(601:1806,2))   

    e3 = delCdt(1309:(end-1),2) - decon(1309:(end-1),2); % look at MSE from 1958-present
    misfit3 = e3'*e3/length(delCdt(1309:(end-1),2));  
    
      C = cov(decon(1309:(end-1),2),yhat2(1309:(end-1),1));
   
  [R,P,RLO,RUP] = corrcoef(yhat2(1297:(end-1),1),decon(1297:(end-1),2));
  
  %R(1,2)^2

  error1= betahat(1)-ci(1);
    
figure
plot(decon(:,1),decon(:,2),delCdt(:,1),yhat2)
xlabel('year')
ylabel('ppm CO2/year')
title('land uptake')
legend('Residual uptake','land uptake with T effects','location','northwest')
grid
end

    


% Do "reverse deconvolution" to calculate modeled atmospheric change in
% CO2
i4 = find(fas(:,1) == start_year);
j4 = find(fas(:,1) == end_year);

%i5 = find(fas(:,1) == start_year);
%j5 = find(fas(:,1) == end_year);

if(LUlevel==1)
    newat(:,1) = year2;
    newat(:,2) =  ff(:,2) - Aoc*fas(i4:j4,2) + LU(:,2) + delCdt(:,2);
    
    % JLD code below:
%     newat(:,2) =  ff(i2:j2,2)....
%     - Aoc*fas(i3:j3,2) + LU(1:length(year2),2) + delCdt(:,2) ;

elseif(LUlevel==2)
    
    % 1850 to 2005  --- units on this?
    newat(:,1) = year2;
newat(:,2) =  ff(:,2) - Aoc*fas(i4:end,2) + LUex(:,2) + delCdt(:,2);
%JLD 3/13/18 changed from delCdt(:,2) to delCdt(1:1867,2)
%JLD 3/21/18 this case will continue past where LR plots it- should match
%up to 2005.5

    % JLD code below:
%     newat(:,1) = year2;
%     newat(:,2) =  ff(i2:j2,2)....
%     - Aoc*fas(i3:j3,2) + LUex(1:length(year2),2) + delCdt(:,2) ;
end

% combine CO2a from mlointerp (goes through 2010) with yearly thru 2016
% starts at year 1 (?) incremented by year
% CO2_2016 = csvread('mergedCO2_2016.csv'); 
% CO2_2016mo(:,1) = year2;
% CO2_2016mo(:,2) = (interp1(CO2_2016(:,1),CO2_2016(:,2),year2)).';
% k1 = find(floor(100*CO2_2016mo(:,1)) == floor(100*CO2a(end,1)));
    


atmcalc(:,1) = year2;
atmcalc(:,2) = cumsum(newat(:,2)/12);
atmcalc(:,2) = atmcalc(:,2)+co2_preind;

co2_diff(:,1) = year2;
co2_diff(:,2) = CO2a(:,2)-atmcalc(:,2);
i6 = find(co2_diff(:,1) == 1959);
j6 = find(co2_diff(:,1) == 1979);
meandiff = mean(co2_diff(i6:j6,2)); % mean difference over 1959-1979
atmcalc2 = atmcalc(:,2)+meandiff;

% co2_diff(:,1) = year2;
% co2_diff(:,2) = CO2_2016mo(:,2)-atmcalc(:,2);
% i6 = find(co2_diff(:,1) == 1959);
% j6 = find(co2_diff(:,1) == 1979);
% meandiff = mean(co2_diff(i6:j6,2)); % mean difference over 1959-1979
% atmcalc2 = atmcalc(:,2)+meandiff;

figure
%plot(year2,atmcalc2, CO2_2016mo(:,1), CO2_2016mo(:,2));
plot(year2,atmcalc2,CO2a(:,1),CO2a(:,2));
xlabel('year')
ylabel('ppm CO2')
title('modeled atmospheric co2 from deconvolution')
legend('modeled atmospheric co2','observed atmos','location','northwest');
grid


elseif strcmpi('no',inputStr2)
    disp('all done!')
end


