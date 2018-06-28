% error_calcs.m
%
% june 15, 2018
%
% julia dohner



%% load data

addpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle');
addpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle/matlab_variables');

clear all

load tempDep;
load tempIndep;

load timeframe2005_ddtvars;
load const2_ddtvars;

%% calculate root mean square error - over optimization period 1900-2015.5

i = find(CconstLU_ddt(:,1) == 1900);

Vhough_RMSEopt = sqrt(mean(Vhough_ddt(i:end,2).^2));
Vhansis_RMSEopt = sqrt(mean(Vhansis_ddt(i:end,2).^2));
VLRLU_RMSEopt = sqrt(mean(VLRLU_ddt(i:end,2).^2));
VLRLUex_RMSEopt = sqrt(mean(VLRLUex_ddt(i:end,2).^2));
Vconst_RMSEopt = sqrt(mean(VconstLU_ddt(i:end,2).^2));
V2005_hough_RMSEopt = sqrt(mean(V2005_hough_ddt(i:end,2).^2));
V2005_hansis_RMSEopt = sqrt(mean(V2005_hansis_ddt(i:end,2).^2));
V2005_hough03_RMSEopt = sqrt(mean(V2005_hough03_ddt(i:end,2).^2));
V2005_const_RMSEopt = sqrt(mean(V2005_const_ddt(i:end,2).^2));
Vpresent_const2_RMSEopt = sqrt(mean(Vpresent_const2_ddt(i:end,2).^2));
V2005_const2_RMSEopt = sqrt(mean(V2005_const2_ddt(i:end,2).^2));

Chough_RMSEopt = sqrt(mean(Chough_ddt(i:end,2).^2));
Chansis_RMSEopt = sqrt(mean(Chansis_ddt(i:end,2).^2));
CLRLU_RMSEopt = sqrt(mean(CLRLU_ddt(i:end,2).^2));
CLRLUex_RMSEopt = sqrt(mean(CLRLUex_ddt(i:end,2).^2));
Cconst_RMSEopt = sqrt(mean(CconstLU_ddt(i:end,2).^2));
C2005_hough_RMSEopt = sqrt(mean(C2005_hough_ddt(i:end,2).^2));
C2005_hansis_RMSEopt = sqrt(mean(C2005_hansis_ddt(i:end,2).^2));
C2005_hough03_RMSEopt = sqrt(mean(C2005_hough03_ddt(i:end,2).^2));
C2005_const_RMSEopt = sqrt(mean(C2005_const_ddt(i:end,2).^2));
Cpresent_const2_RMSEopt = sqrt(mean(Cpresent_const2_ddt(i:end,2).^2));
C2005_const2_RMSEopt = sqrt(mean(C2005_const2_ddt(i:end,2).^2));


% quick Peters comparison stuff
j = find(CconstLU_ddt(:,1) == 1959);
Vhough_RMSEopt = sqrt(mean(Vhough_ddt(j:end,2).^2));
Vhansis_RMSEopt = sqrt(mean(Vhansis_ddt(j:end,2).^2));
VLRLU_RMSEopt = sqrt(mean(VLRLU_ddt(j:end,2).^2));
VLRLUex_RMSEopt = sqrt(mean(VLRLUex_ddt(j:end,2).^2));
Vconst_RMSEopt = sqrt(mean(VconstLU_ddt(j:end,2).^2));
Vconst2_RMSEopt = sqrt(mean(Vpresent_const2_ddt(j:end,2).^2));

meanbud_hough = mean(Vhough_ddt(j:end,2));
meanbud_hansis = mean(Vhansis_ddt(j:end,2));
meanbud_hough03 = mean(VLRLU_ddt(j:end,2));
meanbud_const = mean(VconstLU_ddt(j:end,2));
meanbud_const2 = mean(Vpresent_const2_ddt(j:end,2));

%% 1850 - present

Vhough_RMSEfull = sqrt(mean(Vhough_ddt(:,2).^2));
Vhansis_RMSEfull = sqrt(mean(Vhansis_ddt(:,2).^2));
VLRLU_RMSEfull = sqrt(mean(VLRLU_ddt(:,2).^2));
VLRLUex_RMSEfull = sqrt(mean(VLRLUex_ddt(:,2).^2));
Vconst_RMSEfull = sqrt(mean(VconstLU_ddt(:,2).^2));
V2005_hough_RMSEfull = sqrt(mean(V2005_hough_ddt(:,2).^2));
V2005_hansis_RMSEfull = sqrt(mean(V2005_hansis_ddt(:,2).^2));
V2005_hough03_RMSEfull = sqrt(mean(V2005_hough03_ddt(:,2).^2));
V2005_const_RMSEfull = sqrt(mean(V2005_const_ddt(:,2).^2));
Vpresent_const2_RMSEfull = sqrt(mean(Vpresent_const2_ddt(:,2).^2));
V2005_const2_RMSEfull = sqrt(mean(V2005_const2_ddt(:,2).^2));

Chough_RMSEfull = sqrt(mean(Chough_ddt(:,2).^2));
Chansis_RMSEfull = sqrt(mean(Chansis_ddt(:,2).^2));
CLRLU_RMSEfull = sqrt(mean(CLRLU_ddt(:,2).^2));
CLRLUex_RMSEfull = sqrt(mean(CLRLUex_ddt(:,2).^2));
Cconst_RMSEfull = sqrt(mean(CconstLU_ddt(:,2).^2));
C2005_hough_RMSEfull = sqrt(mean(C2005_hough_ddt(:,2).^2));
C2005_hansis_RMSEfull = sqrt(mean(C2005_hansis_ddt(:,2).^2));
C2005_hough03_RMSEfull = sqrt(mean(C2005_hough03_ddt(:,2).^2));
C2005_const_RMSEfull = sqrt(mean(C2005_const_ddt(:,2).^2));
Cpresent_const2_RMSEfull = sqrt(mean(Cpresent_const2_ddt(:,2).^2));
C2005_const2_RMSEfull = sqrt(mean(C2005_const2_ddt(:,2).^2));


%% unfiltered 

j = find(Vhough_unfiltddt(:,1) == 1900);

Vhough_unfilt_RMSEopt = sqrt(mean(Vhough_unfiltddt(j:end,2).^2));
Vhansis_unfilt_RMSEopt = sqrt(mean(Vhansis_unfiltddt(j:end,2).^2));
VLRLU_unfilt_RMSEopt = sqrt(mean(VLRLU_unfiltddt(j:end,2).^2));
VLRLUex_unfilt_RMSEopt = sqrt(mean(VLRLUex_unfiltddt(j:end,2).^2));
Vconst_unfilt_RMSEopt = sqrt(mean(Vconst_unfiltddt(j:end,2).^2));
V2005_houghUnfilt_RMSEopt = sqrt(mean(V2005_hough_ddtUnfilt(j:end,2).^2));
V2005_hansisUnfilt_RMSEopt = sqrt(mean(V2005_hansis_ddtUnfilt(j:end,2).^2));
V2005_hough03Unfilt_RMSEopt = sqrt(mean(V2005_hough03_ddtUnfilt(j:end,2).^2));
V2005_constUnfilt_RMSEopt = sqrt(mean(V2005_const_ddtUnfilt(j:end,2).^2));
Vpresent_const2Unfilt_RMSEopt = sqrt(mean(Vpresent_const2_ddtUnfilt(j:end,2).^2));
V2005_const2Unfilt_RMSEopt = sqrt(mean(V2005_const2_ddtUnfilt(j:end,2).^2));

Chough_unfilt_RMSEopt = sqrt(mean(Chough_unfiltddt(j:end,2).^2));
Chansis_unfilt_RMSEopt = sqrt(mean(Chansis_unfiltddt(j:end,2).^2));
CLRLU_unfilt_RMSEopt = sqrt(mean(CLRLU_unfiltddt(j:end,2).^2));
CLRLUex_unfilt_RMSEopt = sqrt(mean(CLRLUex_unfiltddt(j:end,2).^2));
Cconst_unfilt_RMSEopt = sqrt(mean(Cconst_unfiltddt(j:end,2).^2));
C2005_houghUnfilt_RMSEopt = sqrt(mean(C2005_hough_ddtUnfilt(j:end,2).^2));
C2005_hansisUnfilt_RMSEopt = sqrt(mean(C2005_hansis_ddtUnfilt(j:end,2).^2));
C2005_hough03Unfilt_RMSEopt = sqrt(mean(C2005_hough03_ddtUnfilt(j:end,2).^2));
C2005_constUnfilt_RMSEopt = sqrt(mean(C2005_const_ddtUnfilt(j:end,2).^2));
Cpresent_const2Unfilt_RMSEopt = sqrt(mean(Cpresent_const2_ddtUnfilt(j:end,2).^2));
C2005_const2Unfilt_RMSEopt = sqrt(mean(C2005_const2_ddtUnfilt(j:end,2).^2));

Vhough_unfilt_RMSEfull = sqrt(mean(Vhough_unfiltddt(:,2).^2));
Vhansis_unfilt_RMSEfull = sqrt(mean(Vhansis_unfiltddt(:,2).^2));
VLRLU_unfilt_RMSEfull = sqrt(mean(VLRLU_unfiltddt(:,2).^2));
VLRLUex_unfilt_RMSEfull = sqrt(mean(VLRLUex_unfiltddt(:,2).^2));
Vconst_unfilt_RMSEfull = sqrt(mean(Vconst_unfiltddt(:,2).^2));
V2005_houghUnfilt_RMSEfull = sqrt(mean(V2005_hough_ddtUnfilt(:,2).^2));
V2005_hansisUnfilt_RMSEfull = sqrt(mean(V2005_hansis_ddtUnfilt(:,2).^2));
V2005_hough03Unfilt_RMSEfull = sqrt(mean(V2005_hough03_ddtUnfilt(:,2).^2));
V2005_constUnfilt_RMSEfull = sqrt(mean(V2005_const_ddtUnfilt(:,2).^2));
Vpresent_const2Unfilt_RMSEfull = sqrt(mean(Vpresent_const2_ddtUnfilt(:,2).^2));
V2005_const2Unfilt_RMSEfull = sqrt(mean(V2005_const2_ddtUnfilt(:,2).^2));


Chough_unfilt_RMSEfull = sqrt(mean(Chough_unfiltddt(:,2).^2));
Chansis_unfilt_RMSEfull = sqrt(mean(Chansis_unfiltddt(:,2).^2));
CLRLU_unfilt_RMSEfull = sqrt(mean(CLRLU_unfiltddt(:,2).^2));
CLRLUex_unfilt_RMSEfull = sqrt(mean(CLRLUex_unfiltddt(:,2).^2));
Cconst_unfilt_RMSEfull = sqrt(mean(Cconst_unfiltddt(:,2).^2));
C2005_houghUnfilt_RMSEfull = sqrt(mean(C2005_hough_ddtUnfilt(:,2).^2));
C2005_hansisUnfilt_RMSEfull = sqrt(mean(C2005_hansis_ddtUnfilt(:,2).^2));
C2005_hough03Unfilt_RMSEfull = sqrt(mean(C2005_hough03_ddtUnfilt(:,2).^2));
C2005_constUnfilt_RMSEfull = sqrt(mean(C2005_const_ddtUnfilt(:,2).^2));
Cpresent_const2Unfilt_RMSEfull = sqrt(mean(Cpresent_const2_ddtUnfilt(:,2).^2));
C2005_const2Unfilt_RMSEfull = sqrt(mean(C2005_const2_ddtUnfilt(:,2).^2));

%% create output tables, save variables

% rmse over optimization period (1900-2015.5)
VrmseOpt_table = {'Vhough','Vhansis','Vhough03','VconstLU',...
    'V2005_hough','V2005_hansis','V2005_hough03','V2005_const',...
    'Vpresent_const2','V2005_const2'; ...
    Vhough_RMSEopt, Vhansis_RMSEopt, VLRLU_RMSEopt, Vconst_RMSEopt, ...
    V2005_hough_RMSEopt,V2005_hansis_RMSEopt,V2005_hough03_RMSEopt,...
    V2005_const_RMSEopt,Vpresent_const2_RMSEopt,V2005_const2_RMSEopt}

CrmseOpt_table = {'Chough','Chansis','Chough03','CconstLU',...
    'C2005_hough','C2005_hansis','C2005_hough03','C2005_const',...
    'Cpresent_const2','C2005_const2'; ...
    Chough_RMSEopt, Chansis_RMSEopt, CLRLU_RMSEopt, Cconst_RMSEopt,...
    C2005_hough_RMSEopt,C2005_hansis_RMSEopt,C2005_hough03_RMSEopt,...
    C2005_const_RMSEopt, Cpresent_const2_RMSEopt,C2005_const2_RMSEopt}

VrmseOpt_unfilt_table = ...
    {'Vhough_unfilt','Vhansis_unfilt','Vhough03_unfilt','Vconst_unfilt',... 
    'V2005_hough_unfilt','V2005_hansis_unfilt','V2005_hough03_unfilt',...
    'V2005_const_unfilt','Vpresent_const2_unfilt','V2005_const2_unfilt';...
    Vhough_unfilt_RMSEopt, Vhansis_unfilt_RMSEopt,...
    VLRLU_unfilt_RMSEopt, Vconst_unfilt_RMSEopt, ...
    V2005_houghUnfilt_RMSEopt, V2005_hansisUnfilt_RMSEopt,...
    V2005_hough03Unfilt_RMSEopt, V2005_constUnfilt_RMSEopt,...
    Vpresent_const2Unfilt_RMSEopt,V2005_const2Unfilt_RMSEopt}

CrmseOpt_unfilt_table = ...
    {'Chough_unfilt','Chansis_unfilt','Chough03_unfilt','Cconst_unfilt',...
    'C2005_hough_unfilt','C2005_hansis_unfilt','C2005_hough03_unfilt',...
    'C2005_const_unfilt','Cpresent_const2_unfilt','C2005_const2_unfilt'; 
    Chough_unfilt_RMSEopt, Chansis_unfilt_RMSEopt,...
    CLRLU_unfilt_RMSEopt, Cconst_unfilt_RMSEopt,...
    C2005_houghUnfilt_RMSEopt, C2005_hansisUnfilt_RMSEopt,...
    C2005_hough03Unfilt_RMSEopt, C2005_constUnfilt_RMSEopt,...
    Cpresent_const2Unfilt_RMSEopt, C2005_const2Unfilt_RMSEopt}


% rmse over full period (1850-2015.5)
VrmseFull_table = {'Vhough','Vhansis','Vhough03','VconstLU',...
    'V2005_hough','V2005_hansis','V2005_hough03','V2005_const',...
    'Vpresent_const2','V2005_const2'; ...
    Vhough_RMSEfull, Vhansis_RMSEfull, VLRLU_RMSEfull, Vconst_RMSEfull, ...
    V2005_hough_RMSEfull,V2005_hansis_RMSEfull,V2005_hough03_RMSEfull,...
    V2005_const_RMSEfull, Vpresent_const2_RMSEfull, V2005_const2_RMSEfull}

 
CrmseFull_table = {'Chough','Chansis','Chough03','CconstLU', ...
    'C2005_hough','C2005_hansis','C2005_hough03','C2005_const',...
    'Cpresent_const2','C2005_const2'; ...
    Chough_RMSEfull, Chansis_RMSEfull, CLRLU_RMSEfull, Cconst_RMSEfull,...
    C2005_hough_RMSEfull,C2005_hansis_RMSEfull,C2005_hough03_RMSEfull,...
    C2005_const_RMSEfull, Cpresent_const2_RMSEfull,C2005_const2_RMSEfull}

 
VrmseFull_unfilt_table = ...
    {'Vhough_unfilt','Vhansis_unfilt','Vhough03_unfilt','Vconst_unfilt',...
        'V2005_hough_unfilt','V2005_hansis_unfilt','V2005_hough03_unfilt',...
        'V2005_const_unfilt','Vpresent_const2_unfilt','V2005_const2_unfilt';...
    Vhough_unfilt_RMSEfull, Vhansis_unfilt_RMSEfull,...
    VLRLU_unfilt_RMSEfull, Vconst_unfilt_RMSEfull, ...
    V2005_houghUnfilt_RMSEfull, V2005_hansisUnfilt_RMSEfull,...
    V2005_hough03Unfilt_RMSEfull, V2005_constUnfilt_RMSEfull,...
    Vpresent_const2Unfilt_RMSEfull,V2005_const2Unfilt_RMSEfull}


CrmseFull_unfilt_table = ...
    {'Chough_unfilt','Chansis_unfilt','Chough03_unfilt','Cconst_unfilt',... 
    'C2005_hough_unfilt','C2005_hansis_unfilt','C2005_hough03_unfilt',...
    'C2005_const_unfilt','Cpresent_const2_unfilt','C2005_const2_unfilt'; 
    Chough_unfilt_RMSEfull, Chansis_unfilt_RMSEfull,...
    CLRLU_unfilt_RMSEfull, Cconst_unfilt_RMSEfull,...
    C2005_houghUnfilt_RMSEfull, C2005_hansisUnfilt_RMSEfull,...
    C2005_hough03Unfilt_RMSEfull, C2005_constUnfilt_RMSEfull,...
    Cpresent_const2Unfilt_RMSEfull,C2005_const2Unfilt_RMSEfull}






save('RMSEfull','VrmseFull_table','CrmseFull_table');
save('RMSEopt','VrmseOpt_table','CrmseOpt_table');

save('RMSEfull_unfilt', 'VrmseFull_unfilt_table','CrmseFull_unfilt_table');

