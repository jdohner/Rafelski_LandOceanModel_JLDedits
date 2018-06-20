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

load ffdata_total;
load LU_records_monthly; % land use datasets at monthly resolution
load constLU_hough;
constLU_hough = LU;

load presentand2005_vars;
load presentand2005_ddt;

% all co2 and residuals for different land use cases
load tempDep; 
load tempIndep;


%% 3-panel (LU)

figure
subplot(3,1,1)
plot(LUhoughmo(:,1),LUhoughmo(:,2),...
    LUhansismo(:,1),LUhansismo(:,2),...
    LUhough03mo(:,1),LUhough03mo(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d);
%     ,...
%     china_totalann(:,1),d2*china_totalann(:,2),'-.',...
%     usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
%     russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
%     WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
%     everythingElse(:,1),d2*everythingElse(:,2),'--',...
%     global_total(:,1),d2*global_total(:,2),':')

%title('Fossil Fuel & Land Use Change Emissions','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton 2017','Hansis 2015','Houghton 2003','Constant'},...
        'Location','Northwest','FontSize', 16, 'NumColumns',2)
%     ,'China','USA','Russia',...
%     'Western Europe & Japan','Total FF excluding W. Europe & Japan',...
%     'Global total FF'},...

xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 6])
yticks([0:6])
grid

subplot(3,1,2);
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhansis_ddt(:,1),Vhansis_ddt(:,2),...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),...
    VconstLU_ddt(:,1), VconstLU_ddt(:,2))
line([1840,2016],[0,0],'linestyle','--');

%title('Variable T: Obs-Modeled CO_2','FontSize', 22);


legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(3,1,3);
plot(Chough_ddt(:,1),Chough_ddt(:,2),...
    Chansis_ddt(:,1),Chansis_ddt(:,2),...
    CLRLU_ddt(:,1),CLRLU_ddt(:,2),...    
    CconstLU_ddt(:,1), CconstLU_ddt(:,2))
line([1840,2016],[0,0],'linestyle','--');
%title('Fixed T: Obs-Modeled CO_2','FontSize', 22);
legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize', 18);
xlabel('year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid


%% 3-panel (FF)

figure
subplot(4,1,[1 2])
plot(china_totalann(:,1),d2*china_totalann(:,2),'-.',...
    usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
    russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
    everythingElse(:,1),d2*everythingElse(:,2),'--',...
    global_total(:,1),d2*global_total(:,2),':')

%title('Fossil Fuel & Land Use Change Emissions','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'China','USA','Russia',...
    'Western Europe & Japan','Total FF excluding W. Europe & Japan',...
    'Global total FF'},...
        'Location','Northwest','FontSize', 16, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 16])
yticks([0:2:16])
grid

subplot(4,1,3);
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhansis_ddt(:,1),Vhansis_ddt(:,2),...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),...
    VconstLU_ddt(:,1), VconstLU_ddt(:,2))
line([1840,2016],[0,0],'linestyle','--');

%title('Variable T: Obs-Modeled CO_2','FontSize', 22);


legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(4,1,4);
plot(Chough_ddt(:,1),Chough_ddt(:,2),...
    Chansis_ddt(:,1),Chansis_ddt(:,2),...
    CLRLU_ddt(:,1),CLRLU_ddt(:,2),...    
    CconstLU_ddt(:,1), CconstLU_ddt(:,2))
line([1840,2016],[0,0],'linestyle','--');
%title('Fixed T: Obs-Modeled CO_2','FontSize', 22);
legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize', 18);
xlabel('year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

%% residuals t dep next to t-indep

figure
subplot(2,1,1)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhansis_ddt(:,1),Vhansis_ddt(:,2),...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),...
    VconstLU_ddt(:,1),VconstLU_ddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Variable T: Obs-Modeled CO_2','FontSize', 22);

legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(2,1,2)
plot(Chough_ddt(:,1),Chough_ddt(:,2),...
    Chansis_ddt(:,1),Chansis_ddt(:,2),...
    CLRLU_ddt(:,1),CLRLU_ddt(:,2),...
    CconstLU_ddt(:,1),CconstLU_ddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Fixed T: Obs-Modeled CO_2','FontSize', 22);

legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize', 18);
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
    LUhough03mo(:,1),LUhough03mo(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d)
%     ,...
%     china_totalann(:,1),d2*china_totalann(:,2),'-.',...
%     usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
%     russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
%     WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
%     everythingElse(:,1),d2*everythingElse(:,2),':',...
%     global_total(:,1),global_total(:,2),'--')

title('Fossil Fuel & Land Use Change Emissions','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
% land use and land cover change emissions in solid lines, ff emissions in
% dash/dot lines, all emissions outs
legend({'Houghton 2017','Hansis 2015','Houghton 2003'},...
        'Location','Northwest','FontSize', 16, 'NumColumns',2)
%     'Constant','China','USA','Russia',...
%     'Western Europe & Japan','Total ff excluding WE & Japan','FF global total'},...

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

legend({'Filtered','Unfiltered'},'location','northwest','FontSize', 18);
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

legend({'Filtered','Unfiltered'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

%% cutting out 10 years of data

load timeframe2005_ddtvars;


figure
subplot(2,2,1)
plot(Vhough_resid(:,1),Vhough_resid(:,2),...
    V2005_hough_resid(:,1),V2005_hough_resid(:,2))
title('Residuals of Different Time Frames - TempDep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton 2017 fit thru 2015.5','Houghton 2017 fit thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-2 8])
grid

subplot(2,2,3)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    V2005_hough_ddt(:,1),V2005_hough_ddt(:,2))
title('Residual Flux of Different Time Frames - TempDep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton 2017 thru 2015.5','Houghton 2017 thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(2,2,2)
plot(Chough_resid(:,1),Chough_resid(:,2),C2005_hough_resid(:,1),C2005_hough_resid(:,2))
title('Residuals of Different Time Frames - TempIndep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton 2017 thru 2015.5','Houghton 2017 thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-2 8])
grid

subplot(2,2,4)
plot(Chough_ddt(:,1),Chough_ddt(:,2),...
    C2005_hough_ddt(:,1),C2005_hough_ddt(:,2))
title('Residual Flux of Different Time Frames - TempIndep','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton 2017 thru 2015.5','Houghton 2017 thru 2005.5'},...
    'Location','Northwest','FontSize', 16);%, 'NumColumns',2)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

%% 3-panel summary plot
figure
subplot(3,1,1)
plot(LUhoughmo(:,1),LUhoughmo(:,2),...
    LUhansismo(:,1),LUhansismo(:,2),...
    LUhough03mo(:,1),LUhough03mo(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d);
%     ,...
%     china_totalann(:,1),d2*china_totalann(:,2),'-.',...
%     usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
%     russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
%     WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
%     everythingElse(:,1),d2*everythingElse(:,2),'--')

title('Fossil Fuel & Land Use Change Emissions','FontSize', 40);
line([1840,2016],[0,0],'linestyle','--');
legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},...
    'Location','Northwest','FontSize', 20, 'NumColumns',2)
%     'China','USA','Russia',...
%     'Western Europe & Japan','All other emissions',},...
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 6])
grid

subplot(3,1,2)
plot(Vhough_ddt(:,1),Vhough_ddt(:,2),...
    Vhansis_ddt(:,1),Vhansis_ddt(:,2),...
    VLRLU_ddt(:,1),VLRLU_ddt(:,2),...
    VconstLU_ddt(:,1),VconstLU_ddt(:,2));
line([1840,2016],[0,0],'linestyle','--');

title('Variable T: Obs-Modeled CO_2','FontSize', 40);

legend({'Houghton 2017','Hansis 2015','Houghton 2003',...
    'Constant'},'location','northwest','FontSize',20);
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

legend({'Filtered','Unfiltered'},'location','northwest','FontSize', 20);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid