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

load presentand2005_vars;
load presentand2005_ddt;

%% calculate root mean square error - over optimization period 1900-2015.5

i = find(CconstLU_ddt(:,1) == 1900);

Vhough_RMSEopt = sqrt(mean(Vhough_ddt(i:end,2).^2));
Vhansis_RMSEopt = sqrt(mean(Vhansis_ddt(i:end,2).^2));
VLRLU_RMSEopt = sqrt(mean(VLRLU_ddt(i:end,2).^2));
VLRLUex_RMSEopt = sqrt(mean(VLRLUex_ddt(i:end,2).^2));
Vconst_RMSEopt = sqrt(mean(VconstLU_ddt(i:end,2).^2));

Chough_RMSEopt = sqrt(mean(Chough_ddt(i:end,2).^2));
Chansis_RMSEopt = sqrt(mean(Chansis_ddt(i:end,2).^2));
CLRLU_RMSEopt = sqrt(mean(CLRLU_ddt(i:end,2).^2));
CLRLUex_RMSEopt = sqrt(mean(CLRLUex_ddt(i:end,2).^2));
Cconst_RMSEopt = sqrt(mean(CconstLU_ddt(i:end,2).^2));

%% 1850 - present

Vhough_RMSEfull = sqrt(mean(Vhough_ddt(:,2).^2));
Vhansis_RMSEfull = sqrt(mean(Vhansis_ddt(:,2).^2));
VLRLU_RMSEfull = sqrt(mean(VLRLU_ddt(:,2).^2));
VLRLUex_RMSEfull = sqrt(mean(VLRLUex_ddt(:,2).^2));
Vconst_RMSEfull = sqrt(mean(VconstLU_ddt(:,2).^2));

Chough_RMSEfull = sqrt(mean(Chough_ddt(:,2).^2));
Chansis_RMSEfull = sqrt(mean(Chansis_ddt(:,2).^2));
CLRLU_RMSEfull = sqrt(mean(CLRLU_ddt(:,2).^2));
CLRLUex_RMSEfull = sqrt(mean(CLRLUex_ddt(:,2).^2));
Cconst_RMSEfull = sqrt(mean(CconstLU_ddt(:,2).^2));

%% unfiltered 

j = find(Vhough_unfiltddt(:,1) == 1900);

Vhough_unfilt_RMSEopt = sqrt(mean(Vhough_unfiltddt(j:end,2).^2));
Vhansis_unfilt_RMSEopt = sqrt(mean(Vhansis_unfiltddt(j:end,2).^2));
VLRLU_unfilt_RMSEopt = sqrt(mean(VLRLU_unfiltddt(j:end,2).^2));
VLRLUex_unfilt_RMSEopt = sqrt(mean(VLRLUex_unfiltddt(j:end,2).^2));
Vconst_unfilt_RMSEopt = sqrt(mean(Vconst_unfiltddt(j:end,2).^2));

Chough_unfilt_RMSEopt = sqrt(mean(Chough_unfiltddt(j:end,2).^2));
Chansis_unfilt_RMSEopt = sqrt(mean(Chansis_unfiltddt(j:end,2).^2));
CLRLU_unfilt_RMSEopt = sqrt(mean(CLRLU_unfiltddt(j:end,2).^2));
CLRLUex_unfilt_RMSEopt = sqrt(mean(CLRLUex_unfiltddt(j:end,2).^2));
Cconst_unfilt_RMSEopt = sqrt(mean(Cconst_unfiltddt(j:end,2).^2));


Vhough_unfilt_RMSEfull = sqrt(mean(Vhough_unfiltddt(:,2).^2));
Vhansis_unfilt_RMSEfull = sqrt(mean(Vhansis_unfiltddt(:,2).^2));
VLRLU_unfilt_RMSEfull = sqrt(mean(VLRLU_unfiltddt(:,2).^2));
VLRLUex_unfilt_RMSEfull = sqrt(mean(VLRLUex_unfiltddt(:,2).^2));
Vconst_unfilt_RMSEfull = sqrt(mean(Vconst_unfiltddt(:,2).^2));

Chough_unfilt_RMSEfull = sqrt(mean(Chough_unfiltddt(:,2).^2));
Chansis_unfilt_RMSEfull = sqrt(mean(Chansis_unfiltddt(:,2).^2));
CLRLU_unfilt_RMSEfull = sqrt(mean(CLRLU_unfiltddt(:,2).^2));
CLRLUex_unfilt_RMSEfull = sqrt(mean(CLRLUex_unfiltddt(:,2).^2));
Cconst_unfilt_RMSEfull = sqrt(mean(Cconst_unfiltddt(:,2).^2));


%% create output tables, save variables

% rmse over optimization period (1900-2015.5)
VrmseOpt_table = {'Vhough','Vhansis','VLR_high','VLR_low','VconstLU'; ...
    Vhough_RMSEopt, Vhansis_RMSEopt, VLRLU_RMSEopt, VLRLUex_RMSEopt, Vconst_RMSEopt}

CrmseOpt_table = {'Chough','Chansis','CLR_high','CLR_low','CconstLU'; ...
    Chough_RMSEopt, Chansis_RMSEopt, CLRLU_RMSEopt, CLRLUex_RMSEopt, Cconst_RMSEopt}

VrmseOpt_unfilt_table = ...
    {'Vhough_unfilt','Vhansis_unfilt','VLRLU_unfilt','VLRLUex_unfilt','Vconst_unfilt'; 
    Vhough_unfilt_RMSEopt, Vhansis_unfilt_RMSEopt,...
    VLRLU_unfilt_RMSEopt, VLRLUex_unfilt_RMSEopt, Vconst_unfilt_RMSEopt}

CrmseOpt_unfilt_table = ...
    {'Chough_unfilt','Chansis_unfilt','CLRLU_unfilt','CLRLUex_unfilt','Cconst_unfilt'; 
    Chough_unfilt_RMSEopt, Chansis_unfilt_RMSEopt,...
    CLRLU_unfilt_RMSEopt, CLRLUex_unfilt_RMSEopt, Cconst_unfilt_RMSEopt}


% rmse over full period (1850-2015.5)
VrmseFull_table = {'Vhough','Vhansis','VLR_high','VLR_low','VconstLU'; ...
    Vhough_RMSEfull, Vhansis_RMSEfull, VLRLU_RMSEfull, VLRLUex_RMSEfull, Vconst_RMSEfull}
 
CrmseFull_table = {'Chough','Chansis','CLR_high','CLR_low','CconstLU'; ...
    Chough_RMSEfull, Chansis_RMSEfull, CLRLU_RMSEfull, CLRLUex_RMSEfull, Cconst_RMSEfull};
 
VrmseFull_unfilt_table = ...
    {'Vhough_unfilt','Vhansis_unfilt','VLRLU_unfilt','VLRLUex_unfilt','Vconst_unfilt'; 
    Vhough_unfilt_RMSEfull, Vhansis_unfilt_RMSEfull,...
    VLRLU_unfilt_RMSEfull, VLRLUex_unfilt_RMSEfull, Vconst_unfilt_RMSEfull}

CrmseFull_unfilt_table = ...
    {'Chough_unfilt','Chansis_unfilt','CLRLU_unfilt','CLRLUex_unfilt','Cconst_unfilt'; 
    Chough_unfilt_RMSEfull, Chansis_unfilt_RMSEfull,...
    CLRLU_unfilt_RMSEfull, CLRLUex_unfilt_RMSEfull,Cconst_unfilt_RMSEfull}




save('RMSEfull','VrmseFull_table','CrmseFull_table');
save('RMSEopt','VrmseOpt_table','CrmseOpt_table');

save('RMSEfull_unfilt', 'VrmseFull_unfilt_table','CrmseFull_unfilt_table');

