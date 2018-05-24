% lu_plots2.m
%
% script for plotting different land use data
%
% author: julia dohner
% may 4, 2018
%
% updated script using only full records (no appended)

clear all

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;
year2 = (start_year:(1/ts):end_year)';
d = 2.31; % ppm to PgC conversion factor (formerly 1/2.31 opp direction)
d1 = 0.001; % teragram to petagram conversion factor


addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match'));

% load Houghton LU data through 2016
% note: Houghton gave me 1850-2015 record, used 2016 Houghton value as
% listed in GCPv1.3
LUhough = csvread('HoughLU_perscomm_2016.csv');
luYear = LUhough(1,1):(1/ts):LUhough(end,1);
LUhoughmo_0 = (interp1(LUhough(:,1),LUhough(:,2),luYear)).';
LUhoughmo(:,1) = luYear;
LUhoughmo(:,2) = LUhoughmo_0*d1; %convert from TgC to PgC

% load Hansis LU data 1849-2016
% record used for GCP, not one from Hasis et. al 2015
LUhansis = csvread('Pongratz2016_GCP_meanPasture_peat.csv'); 
luYear2 = LUhansis(1,1):(1/ts):LUhansis(end,1);
LUhansismo_0 = (interp1(LUhansis(:,1),LUhansis(:,2),luYear2)).';
LUhansismo(:,1) = luYear2;
LUhansismo(:,2) = LUhansismo_0; 

% load GCB averaged historical LU data
% lu data, TgC/yr
LUhist = csvread('GCPv1.3_historicalLU2016.csv');
luYear3 = LUhist(1,1):(1/ts):LUhist(end,1);
LUhistmo_0 = (interp1(LUhist(:,1),LUhist(:,2),luYear3)).';
LUhistmo(:,1) = luYear3;
LUhistmo(:,2) = LUhistmo_0; 

% load LR high & low land use 
load LR_LU_LUex2016.mat; % land use records used in Rafelski 2009
LU(:,2) = LU(:,2)*d; % convert from ppm to PgC
LUex(:,2) = LUex(:,2)*d; % convert from ppm to PgC

% plot all LU cases over one another
figure
subplot(2,2,[1 2])
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2), ...
    LUhistmo(:,1),LUhistmo(:,2),LU(:,1),LU(:,2), LUex(:,1),LUex(:,2))
grid
title('Land Use Datasets', 'FontSize', 24)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 14)
xlabel('year', 'FontSize', 24)
ylabel('GtC/year', 'FontSize', 24)

% load residuals from driver using most recent ff data, all different lu
% cases 5 total - all for CHM-V with constant SST
load obsCO2_record.mat;
CO2a(:,2) = CO2a(:,2)*d;
load updatedYear_vec.mat; % gives year2 vector
load updatedHough_outputvars.mat;
hough_co2 = atmcalc2*d;
hough_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
load updatedhansis_outputvars.mat;
hansis_co2 = atmcalc2*d;
hansis_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
load updatedGCPhist_outputvars.mat;
gcp_co2 = atmcalc2*d;
gcp_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
load updatedLRLU_outputvars.mat;
LRLU_co2 = atmcalc2*d;
LRLU_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
load updatedLRLUex_outputvars.mat;
LRLUex_co2 = atmcalc2*d;
LRLUex_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];

% plot modeled co2 records over one another
subplot(2,2,3)
plot(year2,CO2a(:,2),year2,hough_co2,year2,hansis_co2,year2,gcp_co2,year2,...
    LRLU_co2,year2,LRLUex_co2)
grid
title('Modeled and Observed CO_2 Records', 'FontSize', 24)
legend({'Observed CO2','Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 14)
xlabel('year', 'FontSize', 24)
ylabel('GtC', 'FontSize', 24)

% plot residuals over one another
subplot(2,2,4)
plot(year2,hough_resid(:,2),year2,hansis_resid(:,2),year2,gcp_resid(:,2),...
    year2,LRLU_resid(:,2),year2,LRLUex_resid(:,2))
grid
title('Observed - Modeled CO_2', 'FontSize', 24)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'FontSize', 14)
line([year2(1),year2(end)],[0,0],'linestyle','--');
xlabel('year', 'FontSize', 24)
ylabel('GtC', 'FontSize', 24)


% smoothed residuals

% default smoothing
hough_resid_sm0 = smooth(hough_resid(:,2),59);
hansis_resid_sm0 = smooth(hansis_resid(:,2),59);
gcp_resid_sm0 = smooth(gcp_resid(:,2),59);
LRLU_resid_sm0 = smooth(LRLU_resid(:,2),59);
LRLUex_resid_sm0 = smooth(LRLUex_resid(:,2),59);

% rloess smoothing
hough_resid_rloess0 = smooth(hough_resid(:,1),hough_resid(:,2),0.3,'rloess');
hansis_resid_rloess0 = smooth(hansis_resid(:,1),hansis_resid(:,2),0.3,'rloess');
gcp_resid_rloess0 = smooth(gcp_resid(:,1),gcp_resid(:,2),0.3,'rloess');
LRLU_resid_rloess0 = smooth(LRLU_resid(:,1),LRLU_resid(:,2),0.3,'rloess');
LRLUex_resid_rloess0 = smooth(LRLUex_resid(:,1),LRLUex_resid(:,2),0.3,'rloess');

% smoothed residuals as fluxes - taking annual differences (Jan - Jan)
hough_resid_sm1 = hough_resid_sm0(1:12:end);
hansis_resid_sm1 = hansis_resid_sm0(1:12:end);
gcp_resid_sm1 = gcp_resid_sm0(1:12:end);
LRLU_resid_sm1 = LRLU_resid_sm0(1:12:end);
LRLUex_resid_sm1 = LRLUex_resid_sm0(1:12:end);

hough_resid_rloess1 = hough_resid_rloess0(1:12:end);
hansis_resid_rloess1 = hansis_resid_rloess0(1:12:end);
gcp_resid_rloess1 = gcp_resid_rloess0(1:12:end);
LRLU_resid_rloess1 = LRLU_resid_rloess0(1:12:end);
LRLUex_resid_rloess1 = LRLUex_resid_rloess0(1:12:end);
year_sm = year2(1:12:end);

blankVec(:,1) = year_sm;
blankVec(:,2) = 0;
hough_resid_ddt = blankVec;
hansis_resid_ddt = blankVec;
gcp_resid_ddt = blankVec;
LRLU_resid_ddt = blankVec;
LRLUex_resid_ddt = blankVec;

blankVec2(:,1) = year_sm;
blankVec2(:,2) = 0;
hough_residrloess_ddt = blankVec2;
hansis_residrloess_ddt = blankVec2;
gcp_residrloess_ddt = blankVec2;
LRLU_residrloess_ddt = blankVec2;
LRLUex_residrloess_ddt = blankVec2;

% calculating derivative as Jan value - Jan value from previous year
for i = 1:length(year_sm)-1
    
    hough_resid_ddt(i,2) = hough_resid_sm1(i+1)-hough_resid_sm1(i);
    hansis_resid_ddt(i,2) = hansis_resid_sm1(i+1)-hansis_resid_sm1(i);
    gcp_resid_ddt(i,2) = gcp_resid_sm1(i+1)-gcp_resid_sm1(i);
    LRLU_resid_ddt(i,2) = LRLU_resid_sm1(i+1)-LRLU_resid_sm1(i);
    LRLUex_resid_ddt(i,2) = LRLUex_resid_sm1(i+1)-LRLUex_resid_sm1(i);
    
    
    hough_residrloess_ddt(i,2) = hough_resid_rloess1(i+1)-hough_resid_rloess1(i);
    hansis_residrloess_ddt(i,2) = hansis_resid_rloess1(i+1)-hansis_resid_rloess1(i);
    gcp_residrloess_ddt(i,2) = gcp_resid_rloess1(i+1)-gcp_resid_rloess1(i);
    LRLU_residrloess_ddt(i,2) = LRLU_resid_rloess1(i+1)-LRLU_resid_rloess1(i);
    LRLUex_residrloess_ddt(i,2) = LRLUex_resid_rloess1(i+1)-LRLUex_resid_rloess1(i);

end

figure
subplot(3,2,[1 2])
h1a = plot(year2,hansis_resid_sm0);
hold on;
h1b = plot(year2,hansis_resid_rloess0);
h2a = plot(year2, hough_resid_sm0);
h2b = plot(year2,hough_resid_rloess0);
h3a = plot(year2, gcp_resid_sm0);
h3b = plot(year2,gcp_resid_rloess0);
h4a = plot(year2, LRLU_resid_sm0);
h4b = plot(year2,LRLU_resid_rloess0);
h5a = plot(year2, LRLUex_resid_sm0);
h5b = plot(year2,LRLUex_resid_rloess0);
% plot(year2,hansis_resid_sm0,year2,hansis_resid_rloess0,...
%     year2, hough_resid_sm0,year2,hough_resid_rloess0,...
%     year2, gcp_resid_sm0,year2,gcp_resid_rloess0,...
%     year2, LRLU_resid_sm0,year2,LRLU_resid_rloess0,...
%     year2, LRLUex_resid_sm0,year2,LRLUex_resid_rloess0)
line([year2(1),year2(end)],[0,0],'linestyle','--');
title('Smoothed residuals (obs-calc) of land use cases (default, rloess smoothing)',...
    'FontSize', 18);
legend([h1a h3a h2a  h4a h5a],{'hansis','GCP','houghton','Rafelski high','Rafelski low'},...
    'FontSize', 12);
% legend('hansis smoothed','','houghton smoothed','','gcp smoothed','',...
%     'Rafelski high LU smoothed','','Rafelski low LU smoothed','',...
%     'location','northwest')
grid
xlabel('year','FontSize', 18)
ylabel('GtC','FontSize', 18)
subplot(3,2,3)
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),...
    hough_residrloess_ddt(:,1),hough_residrloess_ddt(:,2))
line([year2(1),year2(end)],[0,0],'linestyle','--');
title('Houghton - Derivative of residuals (default smooth, rloess)',...
    'FontSize', 18);
% legend('Houghton - derivative of defalt smoothed residuals',...
%     'Houghton - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 18)
ylabel('GtC/year','FontSize', 18)
subplot(3,2,4)
plot(hansis_resid_ddt(:,1),hansis_resid_ddt(:,2),...
    hansis_residrloess_ddt(:,1),hansis_residrloess_ddt(:,2))
line([year2(1),year2(end)],[0,0],'linestyle','--');
title('Hansis - Derivative of residuals (default smooth, rloess)',...
    'FontSize', 18);
% legend('Hansis - derivative of defalt smoothed residuals',...
%     'Hansis - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 18)
ylabel('GtC/year','FontSize', 18)
subplot(3,2,5)
plot(gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),...
    gcp_residrloess_ddt(:,1),gcp_residrloess_ddt(:,2))
line([year2(1),year2(end)],[0,0],'linestyle','--');
title('GCP - Derivative of residuals (default smooth, rloess',...
    'FontSize', 18);
% legend('GCP - derivative of defalt smoothed residuals',...
%     'GCP - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 18)
ylabel('GtC/year','FontSize', 18)
subplot(3,2,6)
plot(LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),...
    LRLU_residrloess_ddt(:,1),LRLU_residrloess_ddt(:,2),...
    LRLUex_resid_ddt(:,1),LRLUex_resid_ddt(:,2),...
    LRLUex_residrloess_ddt(:,1),LRLUex_residrloess_ddt(:,2))
line([year2(1),year2(end)],[0,0],'linestyle','--');
title('LR high and low - Derivative of residuals (default smooth, rloess)',...
    'FontSize', 18);
% legend('LR high - derivative of defalt smoothed residuals',...
%     'LR high - derivative of rloess smoothed residuals', ...
%     'LR low - derivative of defalt smoothed residuals',...
%     'LR low - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 18)
ylabel('GtC/year','FontSize', 18)
