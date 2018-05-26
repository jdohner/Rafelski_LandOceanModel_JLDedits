% ff_residual_plots.m
%
% overlay of ff data and residuals from LR model
%
% julia dohner
% may 15, 2018

clear all

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_FF_untangle'));
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/plotting'));


% all values in thousand metric tonnes
% multiply by 10e-6 to get to GtC
% multiply by 0.001 to get to million metric tons
% multiple by 3.667 to get from million metric tons to 
load lu_obsDiff_collection.mat;
load top4emitters_vars.mat;
load updatedResid;
load updatedFF;
load updatedLU;
load lu_records;
load china_total2; %last half year cut off
load LU_resid_sm_derivs.mat;
load LU_resid_rloess_derivs.mat;
d2 = 1e-6;


%openfig('landuse_residual_derivatives.fig');

% All emission estimates are expressed in million metric tons of carbon. To
% convert these estimates to units of carbon dioxide (CO2), simply multiply
% these estimates by 3.667.

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
title('derivative of residuals from model');
legend('houghton','hansis','gcp','Rafelski high',...
    'Rafelski low','location','northwest');
xlabel('year')
ylabel('GtC/yr')
xlim([1840 2020])
ylim([-3 3])
grid

subplot(3,1,3)
plot(LUhoughmo(:,1),LUhoughmo(:,2),LUhansismo(:,1),LUhansismo(:,2),...
    LUhistmo(:,1),LUhistmo(:,2),LU(:,1),LU(:,2),LUex(:,1),LUex(:,2))
line([1840,2020],[0,0],'linestyle','--');
title('land use records');
legend('houghton','hansis','gcp','rafelski high','rafelski low','location',...
    'northwest')
xlabel('year')
ylabel('GtC/yr')
xlim([1840 2020])
ylim([-2 4])
grid
% figure
% plot(china_solid(:,1),china_solid(:,2),usa_solid(:,1),usa_solid(:,2), ...
%     india_solid(:,1),india_solid(:,2),russia_solid(:,1),russia_solid(:,2));
% title('solid (coal) emissions by nation');
% legend('china','usa','india','russia','location','northwest');
% ylabel('thousand metric ton carbon');
% xlabel('year');
% xlim([1900 2020])
% grid

% % normalized quantities
% hough_norm =  [updatedFFLU_resid(:,1) updatedFFLU_resid(:,2)./max(updatedFFLU_resid(:,2))];
% china_tot_n = [china_total(:,1) china_total(:,2)./max(china_total(:,2))];
% usa_tot_n = [usa_total(:,1) usa_total(:,2)./max(usa_total(:,2))];
% india_tot_n = [india_total(:,1) india_total(:,2)./max(india_total(:,2))];
% russia_tot_n = [russia_total(:,1) russia_total(:,2)./max(russia_total(:,2))];
% 
% figure
% plot(hough_norm(:,1),hough_norm(:,2),china_tot_n(:,1),china_tot_n(:,2),usa_tot_n(:,1),...
%     usa_tot_n(:,2),india_tot_n(:,1),india_tot_n(:,2),russia_tot_n(:,1),russia_tot_n(:,2));
% title('normalized national emissions with residual overlaid')
% legend('residual','china','usa','india','russia')



% smoothed residuals - usa, china, india, russia



% % default smoothing
% usa_total_sm0 = smooth(usa_total(:,2),5);
% china_total_sm0 = smooth(china_total(:,2),5);
% india_total_sm0 = smooth(india_total(:,2),5);
% russia_total_sm0 = smooth(russia_total(:,2),5);
% 
% % rloess smoothing
% % usa_total_rloess0 = smooth(usa_total(:,1),usa_total(:,2),0.3,'rloess');
% % china_total_rloess0 = smooth(china_total(:,1),china_total(:,2),0.3,'rloess');
% % india_total_rloess0 = smooth(india_total(:,1),india_total(:,2),0.3,'rloess');
% % russia_total_rloess0 = smooth(russia_total(:,1),russia_total(:,2),0.3,'rloess');
% 
% usa_total_rloess0 = usa_total;%smooth(usa_total(:,1),usa_total(:,2),0.3,'rloess');
% china_total_rloess0 = china_total;%smooth(china_total(:,1),china_total(:,2),0.3,'rloess');
% india_total_rloess0 = india_total;%smooth(india_total(:,1),india_total(:,2),0.3,'rloess');
% russia_total_rloess0 = russia_total;%smooth(russia_total(:,1),russia_total(:,2),0.3,'rloess');
% 
% figure
% title('smoothed ff emissions by nation')
% plot(usa_total(:,1),usa_total_sm0,china_total(:,1),china_total_sm0,...
%     india_total(:,1),india_total_sm0,russia_total(:,1),russia_total_sm0);
% hold on
% plot(usa_total(:,1),usa_total_rloess0,china_total(:,1),china_total_rloess0,...
%     india_total(:,1),india_total_rloess0,russia_total(:,1),...
%     russia_total_rloess0);
% xlim([1850 2020])
% legend('usa','china','india','russia','usa','china','india','russia');
% 
% % derivative of ff emissions (rate of change) - annual differences
% usa_ddt = usa_total;
% for i = 1:length(usa_total)-1
%     usa_ddt(i,2) = usa_total_rloess0(i+1,2)-usa_total_rloess0(i);
% end
% 
% china_ddt = china_total;
% for i = 1:length(china_total)-1
%     china_ddt(i,2) = china_total_rloess0(i+1,2)-china_total_rloess0(i);
% end
% 
% india_ddt = india_total;
% for i = 1:length(india_total)-1
%     india_ddt(i,2) = india_total_rloess0(i+1,2)-india_total_rloess0(i);
% end
% 
% russia_ddt = russia_total;
% for i = 1:length(russia_total)-1
%     russia_ddt(i,2) = russia_total_rloess0(i+1,2)-russia_total_rloess0(i);
% end
% 
% figure
% plot(usa_ddt(:,1),usa_ddt(:,2),china_ddt(:,1),china_ddt(:,2),...
%     india_ddt(:,1),india_ddt(:,2),russia_ddt(:,1),russia_ddt(:,2));
% line([1850,2014],[0,0],'linestyle','--');
% xlim([1985 2020]);
% legend('usa','china','india','russia')
% 
% figure
% plot(usa_total(:,1),usa_total(:,2),china_total(:,1),china_total(:,2),india_total(:,1),india_total(:,2),russia_total(:,1),russia_total(:,2))
% legend('usa','china','india','russia')

% % plot derivatives of residuals
% yearUpdated = updatedFFLU_resid(:,1);
% updatedResid_ddt = diff([eps;updatedFFLU_resid(:,2)])./diff([eps;updatedFFLU_resid(:,1)]);
% updatedResid_smooth = smooth(updatedFFLU_resid(:,2),15);
% updatedResidsm_ddt = diff([eps;updatedResid_smooth])./diff([eps;updatedFFLU_resid(:,1)]);
% 
% 
% figure
% subplot(2,1,1)
% plot(yearUpdated,updatedResid_ddt)
% title('derivate of residual using updated ff and lu');
% 
% subplot(2,1,2)
% plot(yearUpdated,updatedResidsm_ddt)
% title('derivate of smoothed residual using updated ff and lu');