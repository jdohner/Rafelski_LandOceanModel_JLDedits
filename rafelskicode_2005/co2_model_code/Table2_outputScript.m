% table 2 output script
%
% author: Julia Dohner

% constant SST
% CO2 fertilization only

% table 1 parameters: sst // land use // ocean uptake // temp dependence

% C _ _ - _ for all of these cases (Table 1 labeling convention)

%% LR setup

load land_temp.mat % land temperature records
load npp_T.mat % NPP-weighted temperature record
load landwt_T_2011.mat % land temperature anomaly

nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?
filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3; 

clear year start_year end_year ts

ts = 12; % timesteps per year
start_year = 1850;
end_year = 2006;%2009+(7/12); 
year = start_year:(1/ts):end_year;
year2 = year';

beta = [0.5;2]; % initial guesses for model fit

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year); 

% % Extend land use record by making recent emissions equal to last
% % record
% landusemo(1874:1916,1) = year(1874:1916);
% landusemo(1874:1916,2) = landusemo(1873,2);
% % % 


%Extend extratropical emissions by assuming emissions are zero
extratrop_landmo(1802:length(year2),1) = landusemo(1802:length(year2),1);
extratrop_landmo(1802:length(year2),2) = 0;

%% Calculate residual land uptake
% run to 8/2009
% using high land use emissions
% i = length(year2);


i1 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(start_year+(1/24))));
j1 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(end_year+(1/24))));

i2 = find(ff1(:,1) == start_year);
j2 = find(ff1(:,1) == end_year);

i3 = find(fas(:,1) == start_year);
j3 = find(fas(:,1) == end_year);

residual(:,1) = year2;
residual(:,2) = dtdelpCO2a(i1:j1,2) - ff1(i2:j2,2)....
+ Aoc*fas(i3:j3,2) - landusemo(1:length(year2),2);

% using extratropical emissions only
residual2(:,1) = year2;
residual2(:,2) = dtdelpCO2a(i1:j1,2) - ff1(i2:j2,2)....
+ Aoc*fas(i3:j3,2) - extratrop_landmo(1:length(year2),2);

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

% 10 year moving boxcar average of land temperature
% [avg_temp] = l_boxcar(tland4,10,12,1,2483,1,2);
% 
% avg_temp(1:60,2) = avg_temp(61,2);

%%----------------------------------------------------------------------%%
% Pick the temperature record to use
%%----------------------------------------------------------------------%%

%%----------------------------------------------------------------------%%
% 
% ***Use these for various T tests***

 temp_anom(1:6,1) =  avg_temp(601:606,1); %Jan 1850-May 1850
 temp_anom(1:6,2) = landtglob(1,2); %355 instead of 1, 360 instead of 6
  
 temp_anom(7:1916,1) = landtglob(1:1910,1); % Starts at the year 1850.5. 
 temp_anom(7:1916,2) = landtglob(1:1910,2); % 
 
  
X = temp_anom(:,:);

%% cases script below:

landuse = [1 2]; % high, low land use
oceanUptake = [3 4 5]; % high (3), medium (4), low (5) ocean uptake
tempDepen = [7 8]; % temp-independent, temp-dependent

cases = (combvec(landuse, oceanUptake, tempDepen))';

for i = 1:length(cases(:,1))
    
    % land use
    if ismember(1,cases) == 1
        LU = 1; % landuse = high
    else 
        LU = 2; % landuse = low
    end
        
    % ocean uptake
    if ismember(3,cases) == 1 
        oceanUptake = 1; % oceanuptake = high
        fas(:,2) = fas(:,2).*
    elseif ismember(4,cases) == 1
        oceanUptake = 2; % oceanuptake = medium
    else 
        oceanUptake = 3; % oceanuptake = low
    end
        
    % temperature dependence
    if ismember(7, cases) == 1
        tempDep = 0; % temp-independent
    else 
        tempDep = 1; % temp-dependent
    end
    
    % run all of the model code for each case (each iteration of i)
    % but don't make any plots - only plot by running in comand line at end
    % in each iteration, add a line to an output table of epsilon, q10,
    % error
    

%%----------------------------------------------------------------------%%

% calculate a 10-year running boxcar average of the residual land uptake
% don't use 10 year mean before 1957, because data are already smoothed
% (ice core)
if(LU==1) %high land use
%[residual10] = l_boxcar(residual,10,12,1,length(residual),1,2);
% l_boxcar(func,boxlength,dt,starttime,endtime,datecol,numcol)
[residual10a] = l_boxcar(residual,10,12,1225,length(residual),1,2);
residual10(1:1284,:) = residual(1:1284,:);
residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

decon = residual;
elseif(LU ==2) % low land use
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

        [betahat,resid,J] = nlinfit(X,residual10(601:end,2),'Table2_tempDep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513; change 601 to 1297

    elseif(filter == 2) % fit to unfiltered record

        [betahat,resid,J] = nlinfit(X,decon(1057:1513,2),'Table2_tempDep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513. After 1958: 1297; was on 1309:end

    end
else
        if(filter == 1) % fit to 10-year filtered record

        [betahat,resid,J] = nlinfit(X,residual10(601:end,2),'Table2_tempIndep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513; change 601 to 1297

        elseif(filter == 2) % fit to unfiltered record

        [betahat,resid,J] = nlinfit(X,decon(1057:1513,2),'Table2_tempIndep_land_fit_Qs_annotate',beta); %change 601:end to 1081:1513. After 1958: 1297; was on 1309:end

    end
    
end

%% Look at covariances and correlations between model result and calculated land uptake 
[N,P] = size(X);
% 
covariance = inv(J(1:1177,:)'*J(1:1177,:))*sum(resid(1:1177,:).^2)/(N-P) 

[isize,jsize] = size(covariance);
for i=1:isize
    for j=1:jsize
correlation(i,j) = covariance(i,j)/sqrt(covariance(i,i)*covariance(j,j));
    end
end

%% Get uncertainties of best fit values
ci = nlparci(betahat,resid,J);

%% Redefine values of epsilon, gamma and Q1

if(nitrogen == 1) % N-fertilization case
    epsilon = 0;% betahat(2)
    gamma = betahat(1)
    Q1 = betahat(2)
    Q2 = 1;%betahat(3)
else % co2-fertilization case
    epsilon = betahat(1)%0.79;%
    Q1 = betahat(2)%4.91;
    Q2 = 1;%betahat(2)
end

year2 = year';

%% Run the best fit values in the model again to plot
 [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 
 
%%% Nitrogen%%%

% [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN_annotate(epsilon,Q1,Q2,gamma,ff1,ts,year2,dpCO2a,X);
   
    
delCdt(:,2) = -delCdt(:,2);

%% 10 year moving boxcar average of model result
[delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);

%% Quantify how good the fit is using the mean squared error (MSE)       
% Use filtered or unfiltered yhat
%----------------------------------------------------------------%
%
% Filtered
%
%------------------------------------------------

if(filter == 1)

yhat2 = delC10(:,2);

%yhat2 = yhat2 + (0 - yhat2(1));
i5 = find(year2 == 1855);
    
    e = delC10(i5:end,2) - residual10(i5:end,2); % look at MSE for 1855 to 2000
    
    misfit = e'*e/length(delC10(i5:end,2));  
    
    % 1900 to 2005
    i6 = find(year2 == 1900);
    e2 = delC10(i6:end,2) - residual10(i6:end,2); % look at MSE for 1900 to 2000
    misfit2 = e2'*e2/length(delC10(i6:end,2));  


   error1 = betahat(1)-ci(1);
   
%   error2= betahat(2)-ci(2)

%   error3= betahat(3)-ci(3)
   
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
legend('Residual uptake','land uptake without T effects','land uptake with T effects')

end

    


%% Do "reverse deconvolution" to calculate modeled atmospheric change in
%% CO2
if(LU==1)
newat(:,1) = year2;
% 1850 to end year
% i2 = find(ff1(:,1) == start_year);
% j2 = find(ff1(:,1) == end_year);
% 
% i3 = find(fas(:,1) == start_year);
% j3 = find(fas(:,1) == end_year);
newat(:,2) =  ff1(i2:j2,2)....
- Aoc*fas(i3:j3,2) + landusemo(1:length(year2),2) + delCdt(:,2) ;

elseif(LU==2)
newat(:,1) = year2;
% 1850 to 2005  --- units on this?
newat(:,2) =  ff1(i2:j2,2)....
- Aoc*fas(i3:j3,2) + extratrop_landmo(1:length(year2),2) + delCdt(:,2) ;
end


    
    
end
    