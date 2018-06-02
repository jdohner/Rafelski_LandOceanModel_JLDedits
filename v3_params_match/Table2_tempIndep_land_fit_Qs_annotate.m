% Things to check: end year defined correctly? line 12
% beta variables defined correctly? lines 22-32
% correct model chosen? lines 36-43
% correct time period for yhat? line 55

function yhat = Table2_tempIndep_land_fit_Qs_annotate(beta,X)

% Initial conditions - set according to Joos et al

% ts = 12;
% start_year = 1850;
% end_year = 2006;%2009+(7/12);
% year2 = (start_year:(1/ts):end_year)';
load yearinfo.mat

CO2 = 1; % 1 = CO2 fertilization, 0 = N fertilization model

% Get parameters

%[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year);
[dtdelpCO2a,dpCO2a,~,~,CO2a] = getObservedCO2_2(ts,start_year,end_year);


%[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3;




if CO2 == 1
    % For CO2 fertilization model
    epsilon = beta(1);
    Q1 = 1; %beta(2);
    Q2 = 1;%beta(2);
    
    [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 
else 
    % For N fertilization model
    epsilon = 0;
    gamma = beta(1);
    Q1 = beta(2);
    Q2 = 1;%beta(3);
    [fas,ff,LU,LUex] = getSourceSink3(year2,ts);
    [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN(epsilon,Q1,Q2,gamma,ff(601:end,:),ts,year2,dpCO2a,X);
end 


delCdt(:,2) = -delCdt(:,2);

% 10-year moving boxcar average of model
[delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);


% yhat is the term that is compared to the residual flux in nlinfit. 
% Change the index numbers here and in nonlin_land_Qs_annotate (e.g. line
% 158) to fit to a different time period

yhat = delC10(601:end,2);