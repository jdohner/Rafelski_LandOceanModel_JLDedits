% ff_lu_untangle_plots.m
%
% overlay of ff data and residuals from LR model for temp dep runs
%
% julia dohner
% may 15, 2018

clear all

tempDep = 1; % temp dep runs (1) or temp independent runs (0)

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle'));
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/plotting'));


% all values in thousand metric tonnes
% multiply by 10e-6 to get to GtC
% multiply by 0.001 to get to million metric tons
% multiple by 3.667 to get from million metric tons to 

d = 2.31; % ppm to PgC conversion factor (formerly 1/2.31 opp direction)
d1 = 0.001; % teragram to petagram conversion factor
d2 = 1e-6;

load obsCO2_record.mat;
CO2a(:,2) = CO2a(:,2)*d; % convert to PgC

load top4emitters_vars.mat;
load china_total2; %last half year cut off

load LU_records_monthly; % land use datasets at monthly resolution
load LU_resids_co2.mat; % residuals and calculated co2 records for each LU case
%load lu_resids_smoothed; % residuals smoothed using default
%load lu_resids_rloess; % residuals smoothed using rloess

if tempDep == 1
    % temperature-dependent runs
    load lu_tempDep_resids_sm;
    load lu_tempDep_resids_rloess;
    
    load lu_tempDep_residFluxes_sm;
    load lu_tempDep_residFluxes_rloess;
else
    % temperature-independent runs
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
title('Land Use Datasets', 'FontSize', 24)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 14)
xlabel('year', 'FontSize', 24)
ylabel('GtC/year', 'FontSize', 24)

% plot modeled co2 records over one another
subplot(2,2,3)
plot(year,CO2a(:,2),year,hough_co2,year,hansis_co2,year,gcp_co2,year,...
    LRLU_co2,year,LRLUex_co2)
grid
title('Modeled and Observed CO_2 Records', 'FontSize', 24)
legend({'Observed CO2','Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'location','northwest','FontSize', 14)
xlabel('year', 'FontSize', 24)
ylabel('GtC', 'FontSize', 24)

% plot residuals over one another
subplot(2,2,4)
plot(year,hough_resid(:,2),year,hansis_resid(:,2),year,gcp_resid(:,2),...
    year,LRLU_resid(:,2),year,LRLUex_resid(:,2))
grid
title('Observed - Modeled CO_2', 'FontSize', 24)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low'}, 'FontSize', 14)
line([year(1),year(end)],[0,0],'linestyle','--');
xlabel('year', 'FontSize', 24)
ylabel('GtC', 'FontSize', 24)

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
line([year(1),year(end)],[0,0],'linestyle','--');
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
line([year(1),year(end)],[0,0],'linestyle','--');
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
line([year(1),year(end)],[0,0],'linestyle','--');
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
line([year(1),year(end)],[0,0],'linestyle','--');
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

%% figure 3 - overlay emissions and land use


figure
plot(china_total(:,1),d2*china_total(:,2),'.',usa_total(:,1),d2*usa_total(:,2),'.', ...
    india_total(:,1),d2*india_total(:,2),'.',russia_total(:,1),d2*russia_total(:,2),'.');
title('total emissions by nation');

ylabel('GtC/year');
xlabel('year');
xlim([1800 2020])
grid

hold on
% plot(hough_residrloess_ddt(:,1),hough_residrloess_ddt(:,2),hansis_residrloess_ddt(:,1),...
%     hansis_residrloess_ddt(:,2), gcp_residrloess_ddt(:,1),gcp_residrloess_ddt(:,2),...
%     LRLU_residrloess_ddt(:,1),LRLU_residrloess_ddt(:,2),LRLUex_residrloess_ddt(:,1),...
%     LRLUex_residrloess_ddt(:,2));
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),hansis_resid_ddt(:,1),...
    hansis_resid_ddt(:,2), gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),LRLUex_resid_ddt(:,1),...
    LRLUex_resid_ddt(:,2));

legend('china','usa','india','russia','houghton','hansis','Rafelski high',...
    'Rafelski low','location','northwest');


%% figure 4 - 3-panel plot

figure
subplot(3,1,1)
plot(china_total(:,1),d2*china_total(:,2),usa_total(:,1),d2*usa_total(:,2), ...
    india_total(:,1),d2*india_total(:,2),russia_total(:,1),d2*russia_total(:,2));
title('total emissions by nation');
line([1840,2020],[0,0],'linestyle','--');
legend('china','usa','india','russia','location','northwest')
xlabel('year')
ylabel('GtC/yr')
xlim([1840 2020])
ylim([0 6])
grid

subplot(3,1,2)
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),hansis_resid_ddt(:,1),...
    hansis_resid_ddt(:,2), gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),LRLUex_resid_ddt(:,1),...
    LRLUex_resid_ddt(:,2));
line([1840,2020],[0,0],'linestyle','--');
title('model residual fluxes (obs-modeled)');
legend('houghton','hansis','gcp','Rafelski high',...
    'Rafelski low','location','northwest');
xlabel('year')
ylabel('GtC/yr')
xlim([1840 2020])
ylim([-3 3])
grid

subplot(3,1,3)
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2),...
    LUgcpmo(:,1),LUgcpmo(:,2),LU(:,1),LU(:,2),LUex(:,1),LUex(:,2))
line([1840,2020],[0,0],'linestyle','--');
title('land use records');
legend('houghton','hansis','gcp','rafelski high','rafelski low','location',...
    'northwest')
xlabel('year')
ylabel('GtC/yr')
xlim([1840 2020])
ylim([-2 4])
grid

