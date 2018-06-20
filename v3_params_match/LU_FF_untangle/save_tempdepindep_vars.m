% save_tempdepindep_vars.m
%
% june 15, 2018
%
%julia dohner

clear all

addpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle');
addpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle/matlab_variables');


load LU_resids_co2_tempDep.mat; % residuals and calculated co2 records for each LU case
load lu_tempDep_resids_sm;
load lu_tempDep_resids_rloess;
load lu_tempDep_residFluxes_sm; % smoothed residuals as fluxes
load lu_tempDep_residFluxes_rloess;
load lu_tempDep_residFlux_unfilt; % unsmoothed houghton residual as flux

Vhough_co2 = hough_co2;
Vhough_resid = hough_resid;
Vhansis_co2 = hansis_co2;
Vhansis_resid = hansis_resid;
VLRLU_co2 = LRLU_co2;
VLRLU_resid = LRLU_resid;
VLRLUex_co2 = LRLUex_co2;
VLRLUex_resid = LRLUex_resid;
Vconst_co2 = constLU_co2;
Vconst_resid = constLU_resid;

Vhough_ddt = hough_resid_ddt;
Vhansis_ddt = hansis_resid_ddt;
VLRLU_ddt = LRLU_resid_ddt;
VLRLUex_ddt = LRLUex_resid_ddt;
VconstLU_ddt = constLU_resid_ddt;

Vhough_unfiltddt = hough_ddt_unfilt;
Vhansis_unfiltddt = hansis_ddt_unfilt;
VLRLU_unfiltddt = LRLU_ddt_unfilt;
VLRLUex_unfiltddt = LRLUex_ddt_unfilt;
Vconst_unfiltddt = constLU_ddt_unfilt;



load LU_resids_co2_tempIndep.mat; % residuals and calculated co2 records for each LU case
load lu_tempIndep_resids_sm;
load lu_tempIndep_resids_rloess;
load lu_tempIndep_residFluxes_sm; % smoothed residuals as fluxes
load lu_tempIndep_residFluxes_rloess;
load lu_tempIndep_residFlux_unfilt; % unsmoothed houghton residual as flux

Chough_co2 = hough_co2;
Chough_resid = hough_resid;
Chansis_co2 = hansis_co2;
Chansis_resid = hansis_resid;
CLRLU_co2 = LRLU_co2;
CLRLU_resid = LRLU_resid;
CLRLUex_co2 = LRLUex_co2;
CLRLUex_resid = LRLUex_resid;
Cconst_co2 = constLU_co2;
Cconst_resid = constLU_resid;

Chough_ddt = hough_resid_ddt;
Chansis_ddt = hansis_resid_ddt;
CLRLU_ddt = LRLU_resid_ddt;
CLRLUex_ddt = LRLUex_resid_ddt;
CconstLU_ddt = constLU_resid_ddt;

Chough_unfiltddt = hough_ddt_unfilt;
Chansis_unfiltddt = hansis_ddt_unfilt;
CLRLU_unfiltddt = LRLU_ddt_unfilt;
CLRLUex_unfiltddt = LRLUex_ddt_unfilt;
Cconst_unfiltddt = constLU_ddt_unfilt;

save('tempIndep','Chough_co2','Chough_resid','Chough_ddt',...
    'Chansis_co2','Chansis_resid','Chansis_ddt',...
    'CLRLU_co2','CLRLU_resid','CLRLU_ddt',...
    'CLRLUex_co2','CLRLUex_resid','CLRLUex_ddt',...
    'Cconst_co2','Cconst_resid','CconstLU_ddt',...
    'Chough_unfiltddt','Chansis_unfiltddt',...
    'CLRLU_unfiltddt','CLRLUex_unfiltddt','Cconst_unfiltddt')

save('tempDep','Vhough_co2','Vhough_resid','Vhough_ddt',...
    'Vhansis_co2','Vhansis_resid','Vhansis_ddt',...
    'VLRLU_co2','VLRLU_resid','VLRLU_ddt',...
    'VLRLUex_co2','VLRLUex_resid','VLRLUex_ddt',...
    'Vconst_co2','Vconst_resid','VconstLU_ddt',...
    'Vhough_unfiltddt','Vhansis_unfiltddt',...
    'VLRLU_unfiltddt','VLRLUex_unfiltddt','Vconst_unfiltddt')




