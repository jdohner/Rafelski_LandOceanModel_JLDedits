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
load timeframe2005_vars; % load output for houghton run thru present and 2005
load timeframe_const2; % load output for const2 run (doubling the value from const)

% load residuals from driver using most recent ff data, all different lu
% cases 5 total - all for CHM-V with constant SST
load obsCO2_record.mat;
CO2a(:,2) = CO2a(:,2)*d;
load updatedYear_vec.mat; % gives year2 vector

if tempDep == 1
    load hough_co2_resid_Tdep.mat;
    hough_co2 = atmcalc2*d;
    hough_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load hansis_co2_resid_Tdep.mat;
    hansis_co2 = atmcalc2*d;
    hansis_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load GCP_co2_resid_Tdep.mat;
    gcp_co2 = atmcalc2*d;
    gcp_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load LRLU_co2_resid_Tdep.mat;
    LRLU_co2 = atmcalc2*d;
    LRLU_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load LRLUex_co2_resid_Tdep.mat;
    LRLUex_co2 = atmcalc2*d;
    LRLUex_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load constLU_co2_resid_Tdep.mat;
    constLU_co2 = atmcalc2*d;
    constLU_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];

    % save the co2 residuals for each lu record along with modeled co2 record
    % for each (saving second time so all have appropriate variable names)
    save('LU_resids_co2_tempDep','hough_co2','hough_resid','hansis_co2','hansis_resid',...
        'gcp_co2','gcp_resid','LRLU_co2','LRLU_resid','LRLUex_co2','LRLUex_resid',...
        'constLU_co2','constLU_resid')
else
    load hough_co2_resid_Tindep.mat;
    hough_co2 = atmcalc2*d;
    hough_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load hansis_co2_resid_Tindep.mat;
    hansis_co2 = atmcalc2*d;
    hansis_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load GCP_co2_resid_Tindep.mat;
    gcp_co2 = atmcalc2*d;
    gcp_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load LRLU_co2_resid_Tindep.mat;
    LRLU_co2 = atmcalc2*d;
    LRLU_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load LRLUex_co2_resid_Tindep.mat;
    LRLUex_co2 = atmcalc2*d;
    LRLUex_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];
    load constLU_co2_resid_Tindep.mat;
    constLU_co2 = atmcalc2*d;
    constLU_resid = [obsCalcDiff(:,1) obsCalcDiff(:,2)*d];

    % save the co2 residuals for each lu record along with modeled co2 record
    % for each (saving second time so all have appropriate variable names)
    save('LU_resids_co2_tempIndep','hough_co2','hough_resid','hansis_co2','hansis_resid',...
        'gcp_co2','gcp_resid','LRLU_co2','LRLU_resid','LRLUex_co2','LRLUex_resid',...
        'constLU_co2','constLU_resid')

end


% smoothed residuals

% default smoothing
hough_resid_sm0 = smooth(hough_resid(:,2),59);
hansis_resid_sm0 = smooth(hansis_resid(:,2),59);
gcp_resid_sm0 = smooth(gcp_resid(:,2),59);
LRLU_resid_sm0 = smooth(LRLU_resid(:,2),59);
LRLUex_resid_sm0 = smooth(LRLUex_resid(:,2),59);
constLU_resid_sm0 = smooth(constLU_resid(:,2),59);

% Cpresent_resid_sm0 = smooth(Cpresent_resid(:,2),59);
% Vpresent_resid_sm0 = smooth(Vpresent_resid(:,2),59);
C2005_hough_sm0 = smooth(C2005_hough_resid(:,2).*d,59);
V2005_hough_sm0 = smooth(V2005_hough_resid(:,2).*d,59);
C2005_hansis_sm0 = smooth(C2005_hansis_resid(:,2).*d,59);
V2005_hansis_sm0 = smooth(V2005_hansis_resid(:,2).*d,59);
C2005_hough03_sm0 = smooth(C2005_hough03_resid(:,2).*d,59);
V2005_hough03_sm0 = smooth(V2005_hough03_resid(:,2).*d,59);
C2005_const_sm0 = smooth(C2005_constant_resid(:,2).*d,59);
V2005_const_sm0 = smooth(V2005_constant_resid(:,2).*d,59);

Vpresent_const2_sm0 = smooth(Vpresent_constant2_resid(:,2).*d,59);
Cpresent_const2_sm0 = smooth(Cpresent_constant2_resid(:,2).*d,59);
V2005_const2_sm0 = smooth(V2005_constant2_resid(:,2).*d,59);
C2005_const2_sm0 = smooth(C2005_constant2_resid(:,2).*d,59);


if tempDep == 1
    save('lu_tempDep_resids_sm','hough_resid_sm0','hansis_resid_sm0',...
        'gcp_resid_sm0','LRLU_resid_sm0','LRLUex_resid_sm0',...
        'constLU_resid_sm0')
else 
    save('lu_tempIndep_resids_sm','hough_resid_sm0','hansis_resid_sm0',...
    'gcp_resid_sm0','LRLU_resid_sm0','LRLUex_resid_sm0',...
    'constLU_resid_sm0')
end

% rloess smoothing
hough_resid_rloess0 = smooth(hough_resid(:,1),hough_resid(:,2),0.3,'rloess');
hansis_resid_rloess0 = smooth(hansis_resid(:,1),hansis_resid(:,2),0.3,'rloess');
gcp_resid_rloess0 = smooth(gcp_resid(:,1),gcp_resid(:,2),0.3,'rloess');
LRLU_resid_rloess0 = smooth(LRLU_resid(:,1),LRLU_resid(:,2),0.3,'rloess');
LRLUex_resid_rloess0 = smooth(LRLUex_resid(:,1),LRLUex_resid(:,2),0.3,'rloess');
constLU_resid_rloess0 = smooth(constLU_resid(:,1),constLU_resid(:,2),0.3,'rloess');
% leaving out the timeframe houghton runs from rloess smoothing

if tempDep == 1
    save('lu_tempDep_resids_rloess','hough_resid_rloess0','hansis_resid_rloess0',...
    'gcp_resid_rloess0','LRLU_resid_rloess0','LRLUex_resid_rloess0',...
    'constLU_resid_rloess0')
else
    save('lu_tempIndep_resids_rloess','hough_resid_rloess0','hansis_resid_rloess0',...
        'gcp_resid_rloess0','LRLU_resid_rloess0','LRLUex_resid_rloess0',...
        'constLU_resid_rloess0')
end

% smoothed residuals as fluxes - taking annual differences (Jan - Jan)
hough_resid_sm1 = hough_resid_sm0(1:12:end);
hansis_resid_sm1 = hansis_resid_sm0(1:12:end);
gcp_resid_sm1 = gcp_resid_sm0(1:12:end);
LRLU_resid_sm1 = LRLU_resid_sm0(1:12:end);
LRLUex_resid_sm1 = LRLUex_resid_sm0(1:12:end);
constLU_resid_sm1 = constLU_resid_sm0(1:12:end);
% Vpresent_resid_sm1 = Vpresent_resid_sm0(1:12:end);
% Cpresent_resid_sm1 = Cpresent_resid_sm0(1:12:end);
V2005_hough_sm1 = V2005_hough_sm0(1:12:end);
C2005_hough_sm1 = C2005_hough_sm0(1:12:end);
V2005_hansis_sm1 = V2005_hansis_sm0(1:12:end);
C2005_hansis_sm1 = C2005_hansis_sm0(1:12:end);
V2005_hough03_sm1 = V2005_hough03_sm0(1:12:end);
C2005_hough03_sm1 = C2005_hough03_sm0(1:12:end);
V2005_const_sm1 = V2005_const_sm0(1:12:end);
C2005_const_sm1 = C2005_const_sm0(1:12:end);
Vpresent_const2_sm1 = Vpresent_const2_sm0(1:12:end);
Cpresent_const2_sm1 = Cpresent_const2_sm0(1:12:end);
V2005_const2_sm1 = V2005_const2_sm0(1:12:end);
C2005_const2_sm1 = C2005_const2_sm0(1:12:end);


hough_resid_rloess1 = hough_resid_rloess0(1:12:end);
hansis_resid_rloess1 = hansis_resid_rloess0(1:12:end);
gcp_resid_rloess1 = gcp_resid_rloess0(1:12:end);
LRLU_resid_rloess1 = LRLU_resid_rloess0(1:12:end);
LRLUex_resid_rloess1 = LRLUex_resid_rloess0(1:12:end);
constLU_resid_rloess1 = constLU_resid_rloess0(1:12:end);
year_sm = year2(1:12:end);
year_sm_2005 = V2005_hough_resid(1:12:end,1);

blankVec(:,1) = year_sm(1:end-1);
blankVec(:,2) = zeros(length(year_sm)-1,1);
hough_resid_ddt = blankVec;
hansis_resid_ddt = blankVec;
gcp_resid_ddt = blankVec;
LRLU_resid_ddt = blankVec;
LRLUex_resid_ddt = blankVec;
constLU_resid_ddt = blankVec;
V2005_hough_ddt = blankVec;
C2005_hough_ddt = blankVec;
V2005_hansis_ddt = blankVec;
C2005_hansis_ddt = blankVec;
V2005_hough03_ddt = blankVec;
C2005_hough03_ddt = blankVec;
V2005_const_ddt = blankVec;
C2005_const_ddt = blankVec;
Vpresent_const2_ddt = blankVec;
Cpresent_const2_ddt = blankVec;
V2005_const2_ddt = blankVec;
C2005_const2_ddt = blankVec;

blankVec2(:,1) = year_sm;
blankVec2(:,2) = 0;
hough_residrloess_ddt = blankVec2;
hansis_residrloess_ddt = blankVec2;
gcp_residrloess_ddt = blankVec2;
LRLU_residrloess_ddt = blankVec2;
LRLUex_residrloess_ddt = blankVec2;
constLU_residrloess_ddt = blankVec2;

% calculating derivative as Jan value - Jan value from previous year
for i = 1:length(year_sm)-1
    
    hough_resid_ddt(i,2) = hough_resid_sm1(i+1)-hough_resid_sm1(i);
    hansis_resid_ddt(i,2) = hansis_resid_sm1(i+1)-hansis_resid_sm1(i);
    gcp_resid_ddt(i,2) = gcp_resid_sm1(i+1)-gcp_resid_sm1(i);
    LRLU_resid_ddt(i,2) = LRLU_resid_sm1(i+1)-LRLU_resid_sm1(i);
    LRLUex_resid_ddt(i,2) = LRLUex_resid_sm1(i+1)-LRLUex_resid_sm1(i);
    constLU_resid_ddt(i,2) = constLU_resid_sm1(i+1)-constLU_resid_sm1(i);
    %Vpresent_resid_ddt(i,2) = Vpresent_resid_sm1(i+1)-Vpresent_resid_sm1(i);
    %Cpresent_resid_ddt(i,2) = Cpresent_resid_sm1(i+1)-Cpresent_resid_sm1(i);
    V2005_hough_ddt(i,2) = V2005_hough_sm1(i+1)-V2005_hough_sm1(i);
    C2005_hough_ddt(i,2) = C2005_hough_sm1(i+1)-C2005_hough_sm1(i);
    V2005_hansis_ddt(i,2) = V2005_hansis_sm1(i+1)-V2005_hansis_sm1(i);
    C2005_hansis_ddt(i,2) = C2005_hansis_sm1(i+1)-C2005_hansis_sm1(i);
    V2005_hough03_ddt(i,2) = V2005_hough03_sm1(i+1)-V2005_hough03_sm1(i);
    C2005_hough03_ddt(i,2) = C2005_hough03_sm1(i+1)-C2005_hough03_sm1(i);
    V2005_const_ddt(i,2) = V2005_const_sm1(i+1)-V2005_const_sm1(i);
    C2005_const_ddt(i,2) = C2005_const_sm1(i+1)-C2005_const_sm1(i);
    Vpresent_const2_ddt(i,2) = Vpresent_const2_sm1(i+1)-Vpresent_const2_sm1(i);
    Cpresent_const2_ddt(i,2) = Cpresent_const2_sm1(i+1)-Cpresent_const2_sm1(i);
    V2005_const2_ddt(i,2) = V2005_const2_sm1(i+1)-V2005_const2_sm1(i);
    C2005_const2_ddt(i,2) = C2005_const2_sm1(i+1)-C2005_const2_sm1(i);
    
    hough_residrloess_ddt(i,2) = hough_resid_rloess1(i+1)-hough_resid_rloess1(i);
    hansis_residrloess_ddt(i,2) = hansis_resid_rloess1(i+1)-hansis_resid_rloess1(i);
    gcp_residrloess_ddt(i,2) = gcp_resid_rloess1(i+1)-gcp_resid_rloess1(i);
    LRLU_residrloess_ddt(i,2) = LRLU_resid_rloess1(i+1)-LRLU_resid_rloess1(i);
    LRLUex_residrloess_ddt(i,2) = LRLUex_resid_rloess1(i+1)-LRLUex_resid_rloess1(i);
    constLU_residrloess_ddt(i,2) = constLU_resid_rloess1(i+1)-constLU_resid_rloess1(i);
    
end

% below no longer necessary because fit thru 2005 but plot thru 2015
% blankVec4(:,1) = year_sm_2005; 
% blankVec4(:,2) = 0;
% V2005_resid_ddt = blankVec4;
% C2005_resid_ddt = blankVec4;
% % 2005 runs shorter than all the rest, which are 2015.5
% for i = 1:length(V2005_resid_sm1)-1
%     V2005_resid_ddt(i,2) = V2005_resid_sm1(i+1)-V2005_resid_sm1(i);
%     C2005_resid_ddt(i,2) = C2005_resid_sm1(i+1)-C2005_resid_sm1(i);
% end


blankVec3(:,1) = hough_resid(:,1);
blankVec3(:,2) = 0;

hough_ddt_unfilt = blankVec3;
hansis_ddt_unfilt = blankVec3;
LRLU_ddt_unfilt = blankVec3;
LRLUex_ddt_unfilt = blankVec3;
constLU_ddt_unfilt = blankVec3;
V2005_hough_ddtUnfilt = blankVec3;
C2005_hough_ddtUnfilt = blankVec3;
V2005_hansis_ddtUnfilt = blankVec3;
C2005_hansis_ddtUnfilt = blankVec3;
V2005_hough03_ddtUnfilt = blankVec3;
C2005_hough03_ddtUnfilt = blankVec3;
V2005_const_ddtUnfilt = blankVec3;
C2005_const_ddtUnfilt = blankVec3;
Vpresent_const2_ddtUnfilt = blankVec3;
Cpresent_const2_ddtUnfilt = blankVec3;
V2005_const2_ddtUnfilt = blankVec3;
C2005_const2_ddtUnfilt = blankVec3;

for i = 1:length(hough_resid(:,1))-12
    hough_ddt_unfilt(i,2) = hough_resid(i+12,2)-hough_resid(i,2);
    hansis_ddt_unfilt(i,2) = hansis_resid(i+12,2)-hansis_resid(i,2);
    LRLU_ddt_unfilt(i,2) = LRLU_resid(i+12,2)-LRLU_resid(i,2);
    LRLUex_ddt_unfilt(i,2) = LRLUex_resid(i+12,2)-LRLUex_resid(i,2);    
    constLU_ddt_unfilt(i,2) = constLU_resid(i+12,2)-constLU_resid(i,2);    
    V2005_hough_ddtUnfilt(i,2) = (V2005_hough_resid(i+12,2)-V2005_hough_resid(i,2)).*d;
    C2005_hough_ddtUnfilt(i,2) = (C2005_hough_resid(i+12,2)-V2005_hough_resid(i,2)).*d;
    V2005_hansis_ddtUnfilt(i,2) = (V2005_hansis_resid(i+12,2)-V2005_hansis_resid(i,2)).*d;
    C2005_hansis_ddtUnfilt(i,2) = (C2005_hansis_resid(i+12,2)-V2005_hansis_resid(i,2)).*d;
    V2005_hough03_ddtUnfilt(i,2) = (V2005_hough03_resid(i+12,2)-V2005_hough03_resid(i,2)).*d;
    C2005_hough03_ddtUnfilt(i,2) = (C2005_hough03_resid(i+12,2)-V2005_hough03_resid(i,2)).*d;
    V2005_const_ddtUnfilt(i,2) = (V2005_constant_resid(i+12,2)-V2005_constant_resid(i,2)).*d;
    C2005_const_ddtUnfilt(i,2) = (C2005_constant_resid(i+12,2)-V2005_constant_resid(i,2)).*d;
    Vpresent_const2_ddtUnfilt(i,2) = (Vpresent_constant2_resid(i+12,2)-Vpresent_constant2_resid(i,2)).*d;
    Cpresent_const2_ddtUnfilt(i,2) = (Cpresent_constant2_resid(i+12,2)-Cpresent_constant2_resid(i,2)).*d;
    V2005_const2_ddtUnfilt(i,2) = (V2005_constant2_resid(i+12,2)-V2005_constant2_resid(i,2)).*d;
    C2005_const2_ddtUnfilt(i,2) = (C2005_constant2_resid(i+12,2)-C2005_constant2_resid(i,2)).*d;
end

if tempDep == 1
    save('lu_tempDep_residFluxes_sm','hough_resid_ddt','hansis_resid_ddt',...
        'gcp_resid_ddt','LRLU_resid_ddt','LRLUex_resid_ddt',...
        'constLU_resid_ddt')

    save('lu_tempDep_residFluxes_rloess','hough_residrloess_ddt','hansis_residrloess_ddt',...
        'gcp_residrloess_ddt','LRLU_residrloess_ddt','LRLUex_residrloess_ddt',...
        'constLU_residrloess_ddt')
    
    save('lu_tempDep_residFlux_unfilt','hough_ddt_unfilt','hansis_ddt_unfilt',....
        'LRLU_ddt_unfilt','LRLUex_ddt_unfilt','constLU_ddt_unfilt')
    
    
else % tempDep == 0 (tempIndep)
    save('lu_tempIndep_residFluxes_sm','hough_resid_ddt','hansis_resid_ddt',...
        'gcp_resid_ddt','LRLU_resid_ddt','LRLUex_resid_ddt','constLU_resid_ddt')
        

    save('lu_tempIndep_residFluxes_rloess','hough_residrloess_ddt','hansis_residrloess_ddt',...
        'gcp_residrloess_ddt','LRLU_residrloess_ddt','LRLUex_residrloess_ddt',...
        'constLU_residrloess_ddt')
    
    save('lu_tempIndep_residFlux_unfilt','hough_ddt_unfilt','hansis_ddt_unfilt',....
        'LRLU_ddt_unfilt','LRLUex_ddt_unfilt','constLU_ddt_unfilt')
        
        
end

save('timeframe2005_ddtvars',...
        'V2005_hough_ddt','V2005_hansis_ddt','V2005_hough03_ddt','V2005_const_ddt',...
        'V2005_hough_ddtUnfilt','V2005_hansis_ddtUnfilt','V2005_hough03_ddtUnfilt','V2005_const_ddtUnfilt',...
        'C2005_hough_ddt','C2005_hansis_ddt','C2005_hough03_ddt','C2005_const_ddt',...
        'C2005_hough_ddtUnfilt','C2005_hansis_ddtUnfilt','C2005_hough03_ddtUnfilt','C2005_const_ddtUnfilt')

save('const2_ddtvars',...
        'Vpresent_const2_ddt','Vpresent_const2_ddtUnfilt',...
        'Cpresent_const2_ddt','Cpresent_const2_ddtUnfilt',...
        'V2005_const2_ddt','V2005_const2_ddtUnfilt',...
        'C2005_const2_ddt','C2005_const2_ddtUnfilt')
  