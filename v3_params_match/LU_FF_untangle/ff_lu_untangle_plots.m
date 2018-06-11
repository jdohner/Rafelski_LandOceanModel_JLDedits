% ff_lu_untangle_plots.m
%
% overlay of ff data and residuals from LR model for temp dep runs
%
% julia dohner
% may 15, 2018

clear all

tempDep = 1; % temp dep runs (1) or temp independent runs (0)

addpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle/matlab_variables');

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

%load lu_resids_smoothed; % residuals smoothed using default
%load lu_resids_rloess; % residuals smoothed using rloess

if tempDep == 1
    % temperature-dependent runs
    load LU_resids_co2_tempDep.mat; % residuals and calculated co2 records for each LU case
    
    load lu_tempDep_resids_sm;
    load lu_tempDep_resids_rloess;
    
    load lu_tempDep_residFluxes_sm;
    load lu_tempDep_residFluxes_rloess;
else
    % temperature-independent runs
    load LU_resids_co2_tempIndep.mat; % residuals and calculated co2 records for each LU case
    
    load lu_tempIndep_resids_sm;
    load lu_tempIndep_resids_rloess;
    
    load lu_tempIndep_residFluxes_sm; % smoothed residuals as fluxes
    load lu_tempIndep_residFluxes_rloess;
end




%% figure 1 - all land use cases


% plot all LU cases over one another
figure
subplot(2,2,[1 2])
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2), ...
    LUgcpmo(:,1),LUgcpmo(:,2),LU(:,1),LU(:,2), LUex(:,1),LUex(:,2))
grid
title('Land Use Datasets', 'FontSize', 28)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 18)
xlabel('year', 'FontSize', 28)
set(gca,'FontSize',18)
ylabel('GtC/year', 'FontSize', 28)
set(gca,'FontSize',18)

% plot modeled co2 records over one another
subplot(2,2,3)
plot(year,CO2a(:,2),year,hough_co2,year,hansis_co2,year,gcp_co2,year,...
    LRLU_co2,year,LRLUex_co2)
grid
if tempDep == 1
    title('Temp-dependent: Modeled and Observed CO_2 Records', 'FontSize', 28)
else
    title('Temp-independent: Modeled and Observed CO_2 Records', 'FontSize', 28)
end
legend({'Observed CO2','Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 18)
xlabel('year', 'FontSize', 28)
set(gca,'FontSize',18)
ylabel('GtC', 'FontSize', 28)
set(gca,'FontSize',18)

% plot residuals over one another
subplot(2,2,4)
plot(year,hough_resid(:,2),year,hansis_resid(:,2),year,gcp_resid(:,2),...
    year,LRLU_resid(:,2),year,LRLUex_resid(:,2))
grid
if tempDep ==1 
    title('Temp-dependent: Observed - Modeled CO_2', 'FontSize', 28)
else
    title('Temp-independent: Observed - Modeled CO_2', 'FontSize', 28)
end
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'FontSize', 18)
line([year(1),year(end)],[0,0],'linestyle','--');
xlabel('year', 'FontSize', 28)
set(gca,'FontSize',18)
ylabel('GtC', 'FontSize', 28)
set(gca,'FontSize',18)

%% figure 2 - smoothed residuals as fluxes

figure
subplot(3,2,[1 2])
h1a = plot(year,hansis_resid_sm0);
hold on;
h1b = plot(year,hansis_resid_rloess0);
h2a = plot(year, hough_resid_sm0);
h2b = plot(year,hough_resid_rloess0);
h3a = plot(year, gcp_resid_sm0);
h3b = plot(year,gcp_resid_rloess0);
h4a = plot(year, LRLU_resid_sm0);
h4b = plot(year,LRLU_resid_rloess0);
h5a = plot(year, LRLUex_resid_sm0);
h5b = plot(year,LRLUex_resid_rloess0);
% plot(year,hansis_resid_sm0,year,hansis_resid_rloess0,...
%     year, hough_resid_sm0,year,hough_resid_rloess0,...
%     year, gcp_resid_sm0,year,gcp_resid_rloess0,...
%     year, LRLU_resid_sm0,year,LRLU_resid_rloess0,...
%     year, LRLUex_resid_sm0,year,LRLUex_resid_rloess0)
line([year(1),year(end)],[0,0],'linestyle','--');
if tempDep == 1
    title('Temp-dep: Smoothed residuals (obs-calc) of land use cases (default, rloess smoothing)',...
        'FontSize', 22);
else
    title('Temp-indep: Smoothed residuals (obs-calc) of land use cases (default, rloess smoothing)',...
        'FontSize', 22);
end
legend([h1a h3a h2a  h4a h5a],{'hansis','GCP','houghton','Rafelski high','Rafelski low'},...
    'FontSize', 12);
% legend('hansis smoothed','','houghton smoothed','','gcp smoothed','',...
%     'Rafelski high LU smoothed','','Rafelski low LU smoothed','',...
%     'location','northwest')
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 22)
set(gca,'FontSize',18)
subplot(3,2,3)
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),...
    hough_residrloess_ddt(:,1),hough_residrloess_ddt(:,2))
line([year(1),year(end)],[0,0],'linestyle','--');
title('Houghton - Derivative of residuals (default smooth, rloess)',...
    'FontSize', 22);
% legend('Houghton - derivative of defalt smoothed residuals',...
%     'Houghton - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC/year','FontSize', 22)
set(gca,'FontSize',18)
subplot(3,2,4)
plot(hansis_resid_ddt(:,1),hansis_resid_ddt(:,2),...
    hansis_residrloess_ddt(:,1),hansis_residrloess_ddt(:,2))
line([year(1),year(end)],[0,0],'linestyle','--');
title('Hansis - Derivative of residuals (default smooth, rloess)',...
    'FontSize', 22);
% legend('Hansis - derivative of defalt smoothed residuals',...
%     'Hansis - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC/year','FontSize', 22)
set(gca,'FontSize',18)
subplot(3,2,5)
plot(gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),...
    gcp_residrloess_ddt(:,1),gcp_residrloess_ddt(:,2))
line([year(1),year(end)],[0,0],'linestyle','--');
title('GCP - Derivative of residuals (default smooth, rloess',...
    'FontSize', 22);
% legend('GCP - derivative of defalt smoothed residuals',...
%     'GCP - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC/year','FontSize', 22)
set(gca,'FontSize',18)
subplot(3,2,6)
plot(LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),...
    LRLU_residrloess_ddt(:,1),LRLU_residrloess_ddt(:,2),...
    LRLUex_resid_ddt(:,1),LRLUex_resid_ddt(:,2),...
    LRLUex_residrloess_ddt(:,1),LRLUex_residrloess_ddt(:,2))
line([year(1),year(end)],[0,0],'linestyle','--');
title('LR high and low - Derivative of residuals (default smooth, rloess)',...
    'FontSize', 22);
% legend('LR high - derivative of defalt smoothed residuals',...
%     'LR high - derivative of rloess smoothed residuals', ...
%     'LR low - derivative of defalt smoothed residuals',...
%     'LR low - derivative of rloess smoothed residuals','location',...
%     'northwest')
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC/year','FontSize', 22)
set(gca,'FontSize',18)

%% figure 3 - overlay emissions and land use


figure
plot(china_totalann(:,1),d2*china_totalann(:,2),...
    usa_totalann(:,1),d2*usa_totalann(:,2), ...
    russia_totalann(:,1),d2*russia_totalann(:,2),...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),...
    everythingElse(:,1),d2*everythingElse(:,2),'--');
if tempDep == 1
    title('Fossil Fuel Emissions & Temp-Dependent Land Use Case Residuals','FontSize', 28);
else
    title('Fossil Fuel Emissions & Temp-Independent Land Use Case Residuals','FontSize', 28);
end
ylabel('GtC/year','FontSize', 18);
set(gca,'FontSize',18)
xlabel('Year','FontSize', 18);
set(gca,'FontSize',18)
xlim([1800 2020])
grid

hold on
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),'-.',hansis_resid_ddt(:,1),...
    hansis_resid_ddt(:,2),'-.', gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),'-.',...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),'-.',LRLUex_resid_ddt(:,1),...
    LRLUex_resid_ddt(:,2),'-.');
line([year(1),year(end)],[0,0],'linestyle','--');
legend({'China','USA','Russia','Western Europe & Japan',...
    'Everything else','Houghton','Hansis','GCP','Rafelski high',...
    'Rafelski low'},'location','northwest','FontSize', 18);


%% figure 4 - 3-panel plot

figure
subplot(3,1,1)
plot(china_totalann(:,1),d2*china_totalann(:,2),...
    usa_totalann(:,1),d2*usa_totalann(:,2), ...
    russia_totalann(:,1),d2*russia_totalann(:,2),...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),...
    everythingElse(:,1),d2*everythingElse(:,2),'--');
title('Total Fossil Fuel Emissions by Nation','FontSize', 22);
line([1840,2020],[0,0],'linestyle','--');
legend({'China','USA','Russia','Western Europe & Japan',....
    'All other emissions'},'location','northwest','FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2020])
ylim([0 6])
grid

subplot(3,1,2)
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),hansis_resid_ddt(:,1),...
    hansis_resid_ddt(:,2), gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),LRLUex_resid_ddt(:,1),...
    LRLUex_resid_ddt(:,2));
line([1840,2020],[0,0],'linestyle','--');
if tempDep == 1
    title('Temp-dependent: Model Residual Fluxes (Obs-Modeled)','FontSize', 22);
else
    title('Temp-independent: Model Residual Fluxes (Obs-Modeled)','FontSize', 22);
end

legend({'Houghton','Hansis','GCP','Rafelski high',...
    'Rafelski low'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2020])
ylim([-3 3])
grid

subplot(3,1,3)
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2),...
    LUgcpmo(:,1),LUgcpmo(:,2),LU(:,1),LU(:,2),LUex(:,1),LUex(:,2))
line([1840,2020],[0,0],'linestyle','--');
title('Land Use Records','FontSize', 22);
legend({'Houghton','Hansis','GCP','Rafelski high','Rafelski low'},'location',...
    'northwest','FontSize', 18)
xlabel('year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2020])
ylim([-2 4])
grid

%% t-dep vs. t-indep comparison figure

if tempDep == 1
    hough_co2_tempDep = hough_co2;
    hough_resid_tempDep = hough_resid;
    load LU_resids_co2_tempIndep;
    hough_co2_tempIndep = hough_co2;
    hough_resid_tempIndep = hough_resid;
else
    hough_co2_tempIndep = hough_co2;
    hough_resid_tempIndep = hough_resid;
    load LU_resids_co2_tempDep;
    hough_co2_tempDep = hough_co2;
    hough_resid_tempDep = hough_resid;
end

tempRuns_diff = hough_co2_tempDep - hough_co2_tempIndep;

figure
subplot(3,1,1)
plot(year,hough_co2_tempDep,year,hough_co2_tempIndep)
title('Temp-dependent & Temp-independent runs of Houghton Land Use - Modeled CO_2',...
    'FontSize', 22)
legend({'Hough CO_2 T-dependent','Hough CO_2 T-independent'},...
    'location','northwest','FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
grid

subplot(3,1,2)
plot(year,tempRuns_diff)
line([1840,2020],[0,0],'linestyle','--');
title('Difference between T-dep and T-indep modeled CO_2 using Houghton Land Use',...
    'FontSize', 22)
legend({'Temp-dep - Temp-indep Predicted CO_2'},'FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
grid

subplot(3,1,3)
plot(year,hough_resid_tempDep(:,2),year,hough_resid_tempIndep(:,2))
line([1840,2020],[0,0],'linestyle','--');
title('Temp-dep & Temp-indep runs of Houghton Land Use - Observed-Modeled CO_2 residuals',...
    'FontSize', 22)
legend({'Hough residual T-dependent','Hough residual T-independent'},'FontSize', 14)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
grid

%% 3-panel including constant LU

load constLU_hough_outputvars;
constLU_resid = obsCalcDiff;
constLU_co2 = atmcalc2;

load constLU_hough;
constLU_hough = LU_2016mo;

figure
subplot(3,1,1)
plot(china_totalann(:,1),d2*china_totalann(:,2),...
    usa_totalann(:,1),d2*usa_totalann(:,2), ...
    russia_totalann(:,1),d2*russia_totalann(:,2),...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),...
    everythingElse(:,1),d2*everythingElse(:,2),'--');
title('Total Fossil Fuel Emissions by Nation','FontSize', 22);
line([1840,2020],[0,0],'linestyle','--');
legend({'China','USA','Russia','Western Europe & Japan',....
    'All other emissions'},'location','northwest','FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2020])
ylim([0 6])
grid

subplot(3,1,2)
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),hansis_resid_ddt(:,1),...
    hansis_resid_ddt(:,2), ...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),LRLUex_resid_ddt(:,1),...
    LRLUex_resid_ddt(:,2),constLU_resid_ddt(:,1),constLU_resid_ddt(:,2),'-.r');
line([1840,2020],[0,0],'linestyle','--');
if tempDep == 1
    title('Temp-dependent: Model Residual Fluxes (Obs-Modeled)','FontSize', 22);
else
    title('Temp-independent: Model Residual Fluxes (Obs-Modeled)','FontSize', 22);
end

legend({'Houghton','Hansis','Rafelski high',...
    'Rafelski low','Constant LU'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2020])
ylim([-3 3])
grid

subplot(3,1,3)
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2),...
    LU(:,1),LU(:,2),LUex(:,1),LUex(:,2),...
    constLU_hough(:,1),constLU_hough(:,2)*d)
line([1840,2020],[0,0],'linestyle','--');
title('Land Use Records','FontSize', 22);
legend({'Houghton','Hansis','Rafelski high','Rafelski low','Constant LU'},...
    'location',...
    'northwest','FontSize', 18)
xlabel('year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2020])
ylim([-2 4])
grid

