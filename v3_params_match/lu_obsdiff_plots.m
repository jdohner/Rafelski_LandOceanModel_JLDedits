%lu_obsdiff_plots
%
% author: julia dohner
% may 7, 2018

clear all

load 'lu_obsDiff_collection';

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;
year2 = (start_year:(1/ts):end_year)';

figure
plot(hist_GCB_calc(:,1),hist_GCB_calc(:,2),hist_hough_calc(:,1),hist_hough_calc(:,2),LR_LU_calc(:,1),LR_LU_calc(:,2),LR_LUex_calc(:,1),LR_LUex_calc(:,2))
legend('hist + GCB','hist + houghton','Rafelski 2009 high LU','Rafelski 2009 low LU')
grid
line([year2(1),year2(end)],[0,0]);
ylabel('ppm CO2')
xlabel('year')
title('observed co2 deviation from model')