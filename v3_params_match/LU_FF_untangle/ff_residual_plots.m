% ff_residual_plots.m
%
% overlay of ff data and residuals from LR model
%
% julia dohner
% may 15, 2018

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/plotting'));

% will need to update this once get lu records from Houghton, Pongratz
load lu_obsDiff_collection.mat;
load top4emitters_vars.mat;


%hist_hough_scale = [hist_hough_calc(:,1) 10e.*hist_hough_calc(:,2)];

figure
subplot(3,1,1)
plot(hist_hough_calc(:,1),hist_hough_calc(:,2));
line([1800,2020],[0,0],'linestyle','--');
grid

subplot(3,1,2)
plot(china_total(:,1),china_total(:,2),usa_total(:,1),usa_total(:,2), ...
    india_total(:,1),india_total(:,2),russia_total(:,1),russia_total(:,2));
title('total emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
xlim([1800 2020])
grid

subplot(3,1,3)
plot(china_solid(:,1),china_solid(:,2),usa_solid(:,1),usa_solid(:,2), ...
    india_solid(:,1),india_solid(:,2),russia_solid(:,1),russia_solid(:,2));
title('solid (coal) emissions by nation');
legend('china','usa','india','russia','location','northwest');
ylabel('thousand metric ton carbon');
xlabel('year');
xlim([1900 2020])
grid

% normalized quantities
hough_norm =  [hist_hough_calc(:,1) hist_hough_calc(:,2)./max(hist_hough_calc(:,2))];
china_tot_n = [china_total(:,1) china_total(:,2)./max(china_total(:,2))];
usa_tot_n = [usa_total(:,1) usa_total(:,2)./max(usa_total(:,2))];
india_tot_n = [india_total(:,1) india_total(:,2)./max(india_total(:,2))];
russia_tot_n = [russia_total(:,1) russia_total(:,2)./max(russia_total(:,2))];

figure
plot(hough_norm(:,1),hough_norm(:,2),china_tot_n(:,1),china_tot_n(:,2),usa_tot_n(:,1),...
    usa_tot_n(:,2),india_tot_n(:,1),india_tot_n(:,2),russia_tot_n(:,1),russia_tot_n(:,2));

% need to interpolate these
hough_2015 = csvread('Houghton_global_2015.csv');
figure
plot(hough_2015(:,1),hough_2015(:,2));

% derivatives
hough_ddt = diff(hough_2015(:,2))./diff(hough_2015(:,1));
h_smooth = smooth(hough_2015(:,2),15);
h_sm_ddt = diff(h_smooth)./diff(hough_2015(:,1));
%hough_boxcar = l_boxcar(hist_hough_calc,10,12,1,length(hist_hough_calc),1,2);
