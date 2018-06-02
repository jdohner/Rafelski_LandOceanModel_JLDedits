% getLU_smooth_derivs.m
%
% script for getting smoothed LU records and their derivatives
%
% author: julia dohner
% may 4, 2018
%
% updated script using only full records (no appended)

clear all

tempDep = 1; % indicate whether tempDep (1) or tempIndep (0)
d = 2.31; % ppm to PgC conversion factor (formerly 1/2.31 opp direction)
d1 = 0.001; % teragram to petagram conversion factor

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match'));

load LU_records_monthly; % load LU records at monthly resolution


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

% save the co2 residuals for each lu record along with modeled co2 record
% for each (saving second time so all have appropriate variable names)
save('LU_resids_co2','hough_co2','hough_resid','hansis_co2','hansis_resid',...
    'gcp_co2','gcp_resid','LRLU_co2','LRLU_resid','LRLUex_co2','LRLUex_resid')

% smoothed residuals

% default smoothing
hough_resid_sm0 = smooth(hough_resid(:,2),59);
hansis_resid_sm0 = smooth(hansis_resid(:,2),59);
gcp_resid_sm0 = smooth(gcp_resid(:,2),59);
LRLU_resid_sm0 = smooth(LRLU_resid(:,2),59);
LRLUex_resid_sm0 = smooth(LRLUex_resid(:,2),59);

if tempDep == 1
    save('lu_tempDep_resids_sm','hough_resid_sm0','hansis_resid_sm0',...
        'gcp_resid_sm0','LRLU_resid_sm0','LRLUex_resid_sm0')
else 
    save('lu_tempIndep_resids_sm','hough_resid_sm0','hansis_resid_sm0',...
    'gcp_resid_sm0','LRLU_resid_sm0','LRLUex_resid_sm0')
end

% rloess smoothing
hough_resid_rloess0 = smooth(hough_resid(:,1),hough_resid(:,2),0.3,'rloess');
hansis_resid_rloess0 = smooth(hansis_resid(:,1),hansis_resid(:,2),0.3,'rloess');
gcp_resid_rloess0 = smooth(gcp_resid(:,1),gcp_resid(:,2),0.3,'rloess');
LRLU_resid_rloess0 = smooth(LRLU_resid(:,1),LRLU_resid(:,2),0.3,'rloess');
LRLUex_resid_rloess0 = smooth(LRLUex_resid(:,1),LRLUex_resid(:,2),0.3,'rloess');

if tempDep == 1
    save('lu_tempDep_resids_rloess','hough_resid_rloess0','hansis_resid_rloess0',...
    'gcp_resid_rloess0','LRLU_resid_rloess0','LRLUex_resid_rloess0')
else
    save('lu_tempIndep_resids_rloess','hough_resid_rloess0','hansis_resid_rloess0',...
        'gcp_resid_rloess0','LRLU_resid_rloess0','LRLUex_resid_rloess0')
end

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

blankVec(:,1) = year_sm(1:end-1);
blankVec(:,2) = zeros(length(year_sm)-1,1);
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

if tempDep == 1
    save('lu_tempDep_residFluxes_sm','hough_resid_ddt','hansis_resid_ddt',...
        'gcp_resid_ddt','LRLU_resid_ddt','LRLUex_resid_ddt')

    save('lu_tempDep_residFluxes_rloess','hough_residrloess_ddt','hansis_residrloess_ddt',...
        'gcp_residrloess_ddt','LRLU_residrloess_ddt','LRLUex_residrloess_ddt')
    
else % tempDep == 0 (tempIndep)
    save('lu_tempIndep_residFluxes_sm','hough_resid_ddt','hansis_resid_ddt',...
        'gcp_resid_ddt','LRLU_resid_ddt','LRLUex_resid_ddt')

    save('lu_tempIndep_residFluxes_rloess','hough_residrloess_ddt','hansis_residrloess_ddt',...
        'gcp_residrloess_ddt','LRLU_residrloess_ddt','LRLUex_residrloess_ddt')
        
end
  