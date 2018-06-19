% leftoverPlots.m
% 
% script left over / unused plots from ff_lu_untangle_plots.m
% 
% june 15, 2018
% 
% julia dohner


%% figure 1 - all land use cases


% plot all LU cases over one another
figure
subplot(2,2,[1 2])
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2), ...
    LUgcpmo(:,1),LUgcpmo(:,2),LU(:,1),LU(:,2), ...
    constLU_hough(:,1),constLU_hough(:,2))
grid
title('Land Use Datasets', 'FontSize', 28)
legend({'Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low','Constant'}, 'location','northwest','FontSize', 18)
xlabel('year', 'FontSize', 28)
set(gca,'FontSize',18)
ylabel('GtC/year', 'FontSize', 28)
set(gca,'FontSize',18)

% plot modeled co2 records over one another
subplot(2,2,3)
plot(year,CO2a(:,2),year,hough_co2,year,hansis_co2,year,gcp_co2,year,...
    LRLU_co2,year,LRLUex_co2,year,constLU_co2)
grid
if tempDep == 1
    title('Temp-dependent: Modeled and Observed CO_2 Records', 'FontSize', 28)
else
    title('Temp-independent: Modeled and Observed CO_2 Records', 'FontSize', 28)
end
legend({'Observed CO2','Houghton et al. 2015','Hansis et al. 2017',...
    'Global Carbon Project Average','Rafelski 2009 High',...
    'Rafelski 2009 Low','Constant'}, 'location','northwest','FontSize', 18)
xlabel('year', 'FontSize', 28)
set(gca,'FontSize',18)
ylabel('GtC', 'FontSize', 28)
set(gca,'FontSize',18)

% plot residuals over one another
subplot(2,2,4)
plot(year,hough_resid(:,2),year,hansis_resid(:,2),year,gcp_resid(:,2),...
    year,LRLU_resid(:,2),year,LRLUex_resid(:,2),year,constLU_resid(:,2))
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
h6a = plot(year, constLU_resid_sm0);
h6b = plot(year, constLU_resid_rloess0);
line([year(1),year(end)],[0,0],'linestyle','--');
if tempDep == 1
    title('Temp-dep: Smoothed residuals (obs-calc) of land use cases (default, rloess smoothing)',...
        'FontSize', 22);
else
    title('Temp-indep: Smoothed residuals (obs-calc) of land use cases (default, rloess smoothing)',...
        'FontSize', 22);
end
legend([h1a h3a h2a  h4a h5a h6a],{'hansis','GCP','houghton','Rafelski high',...
    'Rafelski low', 'Constant'},'FontSize', 12);
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
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC/year','FontSize', 22)
set(gca,'FontSize',18)

subplot(3,2,5)
plot(constLU_resid_ddt(:,1),constLU_resid_ddt(:,2),...
    constLU_residrloess_ddt(:,1),constLU_residrloess_ddt(:,2))
line([year(1),year(end)],[0,0],'linestyle','--');
title('Constant - Derivative of residuals (default smooth, rloess',...
    'FontSize', 22);
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
grid
xlabel('year','FontSize', 22)
set(gca,'FontSize',18)
ylabel('GtC/year','FontSize', 22)
set(gca,'FontSize',18)

% figure 3 - overlay emissions and land use


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
xlim([1800 2016])
grid

hold on
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2),'-.',hansis_resid_ddt(:,1),...
    hansis_resid_ddt(:,2),'-.', gcp_resid_ddt(:,1),gcp_resid_ddt(:,2),'-.',...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2),'-.',LRLUex_resid_ddt(:,1),...
    LRLUex_resid_ddt(:,2),'-.',constLU_resid_ddt(:,1),constLU_resid_ddt(:,2),'-.');
line([year(1),year(end)],[0,0],'linestyle','--');
legend({'China','USA','Russia','Western Europe & Japan',...
    'Everything else','Houghton','Hansis','GCP','Rafelski high',...
    'Rafelski low','Constant'},'location','northwest','FontSize', 18);


%% figure 4 - 3-panel plot

figure
subplot(3,1,1)
plot(china_totalann(:,1),d2*china_totalann(:,2),'-.',...
    usa_totalann(:,1),d2*usa_totalann(:,2),'-.', ...
    russia_totalann(:,1),d2*russia_totalann(:,2),'-.',...
    WE_japan_total(:,1),d2*WE_japan_total(:,2),'-.',...
    everythingElse(:,1),d2*everythingElse(:,2),'-*');
title('Total Fossil Fuel Emissions by Nation','FontSize', 22);
line([1840,2016],[0,0],'linestyle','--');
legend({'China','USA','Russia','Western Europe & Japan',....
    'All other emissions'},'location','northwest','FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([0 6])
grid

subplot(3,1,2)
plot(hough_resid_ddt(:,1),hough_resid_ddt(:,2), ...
    hansis_resid_ddt(:,1),hansis_resid_ddt(:,2), ...
    LRLU_resid_ddt(:,1),LRLU_resid_ddt(:,2), ...
    LRLUex_resid_ddt(:,1),LRLUex_resid_ddt(:,2),...
    constLU_resid_ddt(:,1),constLU_resid_ddt(:,2));
    
line([1840,2016],[0,0],'linestyle','--');
if tempDep == 1
    title('Temp-dependent: Model Residual Fluxes (Obs-Modeled)','FontSize', 22);
else
    title('Temp-independent: Model Residual Fluxes (Obs-Modeled)','FontSize', 22);
end

legend({'Houghton','Hansis','GCP','Rafelski high',...
    'Rafelski low','Constant'},'location','northwest','FontSize', 18);
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-3 3])
grid

subplot(3,1,3)
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2),...
    LUgcpmo(:,1),LUgcpmo(:,2),LU(:,1),LU(:,2),LUex(:,1),LUex(:,2),...
    constLU_hough(:,1),constLU_hough(:,2))
line([1840,2016],[0,0],'linestyle','--');
title('Land Use Records','FontSize', 22);
legend({'Houghton','Hansis','GCP','Rafelski high','Rafelski low',....
    'Constant'},'location','northwest','FontSize', 18)
xlabel('year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
xlim([1840 2016])
ylim([-2 4])
grid

load timeframe_presentC;
Cpresent_co2 = atmcalc2;
Cpresent_resid = obsCalcDiff;
Cpresent_q10 = Q1;
Cpresent_eps = epsilon;
Cpresent_year = year2;

load timeframe_2005C;
C2005_co2 = atmcalc2;
C2005_resid = obsCalcDiff;
C2005_q10 = Q1;
C2005_eps = epsilon;
C2005_year = year2;

load timeframe_presentV;
Vpresent_co2 = atmcalc2;
Vpresent_resid = obsCalcDiff;
Vpresent_q10 = Q1;
Vpresent_eps = epsilon;
Vpresent_year = year2;

load timeframe_2005V;
V2005_co2 = atmcalc2;
V2005_resid = obsCalcDiff;
V2005_q10 = Q1;
V2005_eps = epsilon;
V2005_year = year2;

save('presentand2005_vars','Cpresent_co2','Cpresent_resid',...
    'Cpresent_q10','Cpresent_eps','Cpresent_year','C2005_co2',...
    'C2005_resid','C2005_q10','C2005_eps','C2005_year','Vpresent_co2',...
    'Vpresent_resid','Vpresent_q10','Vpresent_eps','Vpresent_year',...
    'V2005_co2','V2005_resid','V2005_q10','V2005_eps','V2005_year')

%% t-dep vs. t-indep comparison figure


tempRuns_diff = hough_ddt_tempDep(:,2) - hough_ddt_tempIndep(:,2);

figure
subplot(3,1,1)
plot(year,hough_co2_tempDep,year,hough_co2_tempIndep, CO2a(:,1),CO2a(:,2),...
    '--')
title('Variable & fixed temp modeled CO_2 - Houghton',...
    'FontSize', 22)
legend({'Variable T','Fixed T','Observed CO_2'},...
    'location','northwest','FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
grid



subplot(3,1,2)
plot(hough_ddt_tempDep(:,1),hough_ddt_tempDep(:,2),hough_ddt_tempIndep(:,1),hough_ddt_tempIndep(:,2))
line([1840,2016],[0,0],'linestyle','--');
title('Temp-dep & Temp-indep residual (obs - modeled) - Houghton',...
    'FontSize', 22)
legend({'Hough residual T-dependent','Hough residual T-independent'},'FontSize', 14)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC/yr','FontSize', 18)
set(gca,'FontSize',18)
grid


subplot(3,1,3)
plot(hough_ddt_tempIndep(:,1),tempRuns_diff)
line([1840,2016],[0,0],'linestyle','--');
title('Difference between T-dep and T-indep residual fluxes - Houghton',...
    'FontSize', 22)
legend({'Temp-dep - Temp-indep Predicted CO_2'},'FontSize', 18)
xlabel('Year','FontSize', 18)
set(gca,'FontSize',18)
ylabel('GtC','FontSize', 18)
set(gca,'FontSize',18)
grid