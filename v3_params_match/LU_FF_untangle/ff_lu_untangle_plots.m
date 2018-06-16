% ff_lu_untangle_plots.m
%
% overlay of ff data and residuals from LR model for temp dep runs
%
% julia dohner
% may 15, 2018

clear all

addpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle/matlab_variables');
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/MATLAB_toolkits'))


% all values in thousand metric tonnes
% multiply by 10e-6 to get to GtC
% multiply by 0.001 to get to million metric tons
% multiple by 3.667 to get from million metric tons to 

d = 2.31; % ppm to PgC conversion factor (formerly 1/2.31 opp direction)
d1 = 0.001; % teragram to petagram conversion factor
d2 = 1e-6;

load obsCO2_record.mat;
CO2a(:,2) = CO2a(:,2)*d; % convert to PgC

%load top4emitters_vars.mat;
%load china_total2; %last half year cut off
load ffdata_total; % 
load LU_records_monthly; % land use datasets at monthly resolution
load constLU_hough;
constLU_hough = LU;

load presentand2005_vars;
load presentand2005_ddt;

% all co2 and residuals for different land use cases
load tempDep; 
load tempIndep;


%% 3-panel including constant LU

figure
%subplot(3,1,1)
h(1) = subplot(3,1,1);


plot(LUhoughmo(:,1),LUhoughmo(:,2),...
    LUhansismo(:,1),LUhansismo(:,2),...
    LU(:,1),LU(:,2),LUex(:,1),LUex(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d,...
    china_totalann(:,1),d2*china_totalann(:,2),'-.',...
    usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
    russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
    everythingElse(:,1),d2*everythingElse(:,2),'--',...
    global_total(:,1),d2*global_total(:,2),':')

%title('Fossil Fuel & Land Use Change Emissions','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant','China','USA','Russia',...
    'Western Europe & Japan','Total FF excluding W. Europe & Japan',...
    'Global total FF'},...
    'Location','Northwest','FontSize', 16, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 6])
yticks([0:6])
grid

h(2) = subplot(3,1,2);
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhansis_ddt(:,1),Vhansis_ddt(:,2),...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),...
    VLRLUex_ddt(:,1),VLRLUex_ddt(:,2),...
    VconstLU_ddt(:,1), VconstLU_ddt(:,2))
line([1840,2016],[0,0],'linestyle','--');

%title('Variable T: Obs-Modeled CO_2','FontSize', 22);


legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant LU'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

h(3) = subplot(3,1,3);
plot(Chough_ddt(:,1),Chough_ddt(:,2),...
    Chansis_ddt(:,1),Chansis_ddt(:,2),...
    CLRLU_ddt(:,1),CLRLU_ddt(:,2),...
    CLRLUex_ddt(:,1),CLRLUex_ddt(:,2),...
    CconstLU_ddt(:,1), CconstLU_ddt(:,2))
line([1840,2016],[0,0],'linestyle','--');
%title('Fixed T: Obs-Modeled CO_2','FontSize', 22);
legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant LU'},'location','northwest','FontSize', 18);
xlabel('year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

% pos1 = get(h(1),'position');
% pos2 = get(h(2),'position');
% pos3 = get(h(3),'position');
% set(h(1), 'position', [l, b, w, h] );



%% residuals t dep next to t-indep

figure
subplot(2,1,1)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),Vhansis_ddt(:,1),...
    Vhansis_ddt(:,2), ...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),VLRLUex_ddt(:,1),...
    VLRLUex_ddt(:,2),VconstLU_ddt(:,1),VconstLU_ddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Variable T: Obs-Modeled CO_2','FontSize', 22);

legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant LU'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(2,1,2)
plot(Chough_ddt(:,1),Chough_ddt(:,2),Chansis_ddt(:,1),...
    Chansis_ddt(:,2), ...
    CLRLU_ddt(:,1),CLRLU_ddt(:,2),CLRLUex_ddt(:,1),...
    CLRLUex_ddt(:,2),CconstLU_ddt(:,1),CconstLU_ddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Fixed T: Obs-Modeled CO_2','FontSize', 22);

legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant LU'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

%% smoothed and unsmoothed derivative flux

figure
subplot(3,1,1)
plot(LUhoughmo(:,1),LUhoughmo(:,2),...
    LUhansismo(:,1),LUhansismo(:,2),...
    LU(:,1),LU(:,2),LUex(:,1),LUex(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d,...
    china_totalann(:,1),d2*china_totalann(:,2),'-.',...
    usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
    russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
    everythingElse(:,1),d2*everythingElse(:,2),':',...
    global_total(:,1),global_total(:,2),'--')

title('Fossil Fuel & Land Use Change Emissions','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
% land use and land cover change emissions in solid lines, ff emissions in
% dash/dot lines, all emissions outs
legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant','China','USA','Russia',...
    'Western Europe & Japan','Total ff excluding WE & Japan','FF global total'},...
    'Location','Northwest','FontSize', 16, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 6])
grid

subplot(3,1,2)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhough_unfiltddt(:,1),Vhough_unfiltddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Variable T: Obs-Modeled CO_2 - Houghton','FontSize', 22);

legend({'Smoothed','Unsmoothed'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(3,1,3)
plot(Chough_ddt(:,1),Chough_ddt(:,2),...
    Chough_unfiltddt(:,1),Chough_unfiltddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Fixed T: Obs-Modeled CO_2 - Houghton','FontSize', 22);

legend({'Smoothed','Unsmoothed'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

%% cutting out 10 years of data




figure
subplot(2,2,1)
plot(Vpresent_year,Vpresent_resid(:,2),V2005_year,V2005_resid(:,2))
title('Residuals of Different Time Frames - TempDep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton thru 2015.5','Houghton thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-2 8])
grid

subplot(2,2,3)
plot(Vpresent_resid_ddt(:,1),Vpresent_resid_ddt(:,2),...
    V2005_resid_ddt(:,1),V2005_resid_ddt(:,2))
title('Residual Flux of Different Time Frames - TempDep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton thru 2015.5','Houghton thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-1 1])
grid

subplot(2,2,2)
plot(Cpresent_year,Cpresent_resid(:,2),C2005_year,C2005_resid(:,2))
title('Residuals of Different Time Frames - TempIndep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton thru 2015.5','Houghton thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-2 8])
grid

subplot(2,2,4)
plot(Cpresent_resid_ddt(:,1),Cpresent_resid_ddt(:,2),...
    V2005_resid_ddt(:,1),V2005_resid_ddt(:,2))
title('Residual Flux of Different Time Frames - TempIndep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton thru 2015.5','Houghton thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-1 1])
grid

%% 3-panel summary plot
figure
subplot(3,1,1)
plot(LUhoughmo(:,1),LUhoughmo(:,2),...
    LUhansismo(:,1),LUhansismo(:,2),...
    LU(:,1),LU(:,2),LUex(:,1),LUex(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d,...
    china_totalann(:,1),d2*china_totalann(:,2),'-.',...
    usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
    russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
    everythingElse(:,1),d2*everythingElse(:,2),'--')

title('Fossil Fuel & Land Use Change Emissions','FontSize', 40);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant','China','USA','Russia',...
    'Western Europe & Japan','All other emissions',},...
    'Location','Northwest','FontSize', 20, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 6])
grid

subplot(3,1,2)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),Vhansis_ddt(:,1),...
    Vhansis_ddt(:,2), ...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),VLRLUex_ddt(:,1),...
    VLRLUex_ddt(:,2),VconstLU_ddt(:,1),VconstLU_ddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Variable T: Obs-Modeled CO_2','FontSize', 40);

legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant LU'},'location','northwest','FontSize',20);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid


subplot(3,1,3)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhough_unfiltddt(:,1),Vhough_unfiltddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Variable T: Obs-Modeled CO_2 - Houghton','FontSize', 40);

legend({'Smoothed','Unsmoothed'},'location','northwest','FontSize', 20);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid