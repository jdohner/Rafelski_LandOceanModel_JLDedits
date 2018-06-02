% lu_plots.m
%
% script for plotting different land use data
%
% author: julia dohner
% may 4, 2018

clear all

ts = 12; % timesteps per year
dt = 1/ts;
start_year = 1850;
end_year = 2015.5;%2009+(7/12);%2015.5;%2011.5;%2009+(7/12); %2006;2015.5; %
year2 = (start_year:(1/ts):end_year)';
d = 1/2.31; % gigaton to ppm conversion factor
month_2016 = (1959:(1/12):2016)';

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/LU_data'));
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/necessary_data'));


% historical LU 1800-1960
hist = csvread('Historical_LU.csv');

% interpolate to monthly
histStart = hist(1,1);
histEnd = hist(end,1);
LU_histmo(:,1) = (histStart:(1/12):histEnd)';
LU_histmo(:,2) = (interp1(hist(:,1),hist(:,2),LU_histmo(:,1)))*d;
histStart = find(LU_histmo(:,1) == year2(1)); % start at 1850
hist = LU_histmo(histStart:end,:); % end at 1960

% LR land use
LU_2006 = csvread('landUse_1800-2006.csv');
LUex_2000 = csvread('landUseExtra_1800-2000.csv');

LU_2006mo(:,1) = (1800:(1/12):2006)';
LU_2006mo(:,2) = interp1(LU_2006(:,1),LU_2006(:,12),LU_2006mo(:,1));
LU_start = find(LU_2006mo(:,1) == year2(1)); % start at 1850
LU_2006mo = LU_2006mo(LU_start:end,:); % end at end year

LUex_2000mo(:,1) = (1800:(1/12):2000)';
LUex_2000mo(:,2) = interp1(LUex_2000(:,1),LUex_2000(:,2),LUex_2000mo(:,1));
LUex_start = find(LUex_2000mo(:,1) == year2(1)); % start at 1850
LUex_2000mo = LUex_2000mo(LUex_start:end,:); % end at end year


% fid = fopen('landUse_1959-2016.txt');
% C = textscan(fid,'%f %f', 'delimiter','\t');
% fclose(fid);
% LU_2016(:,1) = C{1};
% LU_2016(:,2) = C{2};
% % land use 2016 to monthly
% LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),month_2016)).';
% LU_2016mo(:,1) = month_2016;
% LU_2016mo(:,2) = LU_2016mo_0*d;

figure
plot(hist(:,1),hist(:,2));
ylabel('ppm/year')
xlabel('year')

%GCB_LU.csv
%BLUE_LU.csv
%CABLE_LU.csv
%CLASS-CTEM_LU.csv
% global carbon project
GCB = csvread('GCB_LU.csv');
Hough = csvread('Houghton_global_2015.csv');
BLUE = csvread('BLUE_LU.csv');
CABLE = csvread('CABLE_LU.csv');
CLASS = csvread('CLASS-CTEM_LU.csv');


GCBmo(:,1) = month_2016;
GCBmo(:,2) = ((interp1(GCB(:,1),GCB(:,2),month_2016))).*d;
Houghmo(:,1) = month_2016;
Houghmo(:,2) = ((interp1(Hough(:,1),Hough(:,2),month_2016))).*d;
BLUEmo(:,1) = month_2016;
BLUEmo(:,2) = ((interp1(BLUE(:,1),BLUE(:,2),month_2016))).*d;
CABLEmo(:,1) = month_2016;
CABLEmo(:,2) = ((interp1(CABLE(:,1),CABLE(:,2),month_2016))).*d;
CLASSmo(:,1) = month_2016;
CLASSmo(:,2) = ((interp1(CLASS(:,1),CLASS(:,2),month_2016))).*d;

hold on
plot(GCBmo(:,1),GCBmo(:,2),BLUEmo(:,1),BLUEmo(:,2),':',CABLEmo(:,1), ...
    CABLEmo(:,2),':',CLASSmo(:,1),CLASSmo(:,2),':',Houghmo(:,1),Houghmo(:,2)...
    ,LU_2006mo(:,1),LU_2006mo(:,2),'--',LUex_2000mo(:,1),LUex_2000mo(:,2),'--');
legend('historical LU','Global carbon project','BLUE','CABLE','CLASS-CTEM'...
    ,'Houghton','Rafelski 2009','Rafelski extratrop','location','northwest')