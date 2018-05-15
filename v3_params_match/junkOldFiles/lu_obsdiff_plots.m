%lu_obsdiff_plots
%
% author: julia dohner
% may 7, 2018

clear all

load 'lu_obsDiff_collection';
load 'varSST_obsCalcDiff_collection';

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;
year2 = (start_year:(1/ts):end_year)';


figure
plot(hist_GCB_calc(:,1),hist_GCB_calc(:,2),hist_hough_calc(:,1),hist_hough_calc(:,2),LR_LU_calc(:,1),LR_LU_calc(:,2),LR_LUex_calc(:,1),LR_LUex_calc(:,2))
line([year2(1),year2(end)],[0,0]);
legend('hist + GCB','hist + houghton','Rafelski 2009 high LU','Rafelski 2009 low LU')
grid
ylabel('ppm CO2')
xlabel('year')
title('fixed sst - observed co2 deviation from model')

figure
plot(varSST_GCP(:,1),varSST_GCP(:,2),varSST_Hough(:,1),varSST_Hough(:,2),varSST_highLU(:,1),varSST_highLU(:,2),varSST_lowLU(:,1),varSST_lowLU(:,2))
line([year2(1),year2(end)],[0,0]);
legend('hist + GCB','hist + houghton','Rafelski 2009 high LU','Rafelski 2009 low LU')
grid
ylabel('ppm CO2')
xlabel('year')
title('variable sst - observed co2 deviation from model')

figure
plot(varSST_highLU(:,1),varSST_highLU(:,2),LR_LU_calc(:,1),LR_LU_calc(:,2));
line([year2(1),year2(end)],[0,0]);
legend('variable SST','fixed SST')
grid
ylabel('ppm CO2')
xlabel('year')
title('Rafelski 2009 Tropical LU Case - Fixed vs. Variable SST')

figure
plot(varSST_Hough(:,1),varSST_Hough(:,2),hist_hough_calc(:,1),hist_hough_calc(:,2));
line([year2(1),year2(end)],[0,0]);
legend('variable SST','fixed SST')
grid
ylabel('ppm CO2')
xlabel('year')
title('Historical + Houghton LU Case - Fixed vs. Variable SST')