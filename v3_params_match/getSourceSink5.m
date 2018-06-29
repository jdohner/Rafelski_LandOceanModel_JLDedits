% getSourceSink5.m
%
% author: Lauren Rafelski, modified by Julia Dohner
% May 16, 2018
%
% gss version that uses most recent values for FF (Boden 2015, GCP) and 
% land use (Houghton 2015, personal comm.)

function [ff,LU] = getSourceSink5(year, ts,landusedata);

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/v3_params_match/necessary_data'));

d = 1/2.31; % gigaton to ppm conversion factor
d1 = 0.001; % teragram to petagram conversion factor
    
% ff data 1750-2016, GtC/yr
FF_2016 = csvread('GCPv1.3_FF2016.csv');
ffYear = FF_2016(1,1):(1/ts):FF_2016(end,1);
FF_2016mo_0 = (interp1(FF_2016(:,1),FF_2016(:,2),ffYear)).';
FF_2016mo(:,1) = ffYear;
FF_2016mo(:,2) = FF_2016mo_0*d; %convert to pp

if strcmp(landusedata,'hough') || strcmp(landusedata,'const') || strcmp(landusedata,'const2')
    %if using Houghton:
    LU_2016 = csvread('HoughLU_perscomm_2016.csv'); 
    luYear = LU_2016(1,1):(1/ts):LU_2016(end,1);
    LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),luYear)).';
    LU_2016mo(:,1) = luYear;
    %TgC to ppm (d1) only for Houghton
    LU_2016mo(:,2) = LU_2016mo_0*d*d1; %convert from TgC to PgC to ppm
    
    if strcmp(landusedata,'const')
        const = [LU_2016mo(1,1) LU_2016mo(1,2) ; LU_2016mo(end,1) LU_2016mo(end,2)];
        LU_2016mo(:,2) = (interp1(const(:,1),const(:,2),luYear)).';
        i = find(LU_2016mo(:,1) == 1920);
        LU_2016mo(i:end,2) = LU_2016mo(i,2);
    elseif strcmp(landusedata,'const2')
        const2 = [LU_2016mo(1,1) LU_2016mo(1,2) ; LU_2016mo(end,1) 2*LU_2016mo(end,2)];
        LU_2016mo(:,2) = (interp1(const2(:,1),const2(:,2),luYear)).';
        i = find(LU_2016mo(:,1) == 1920);
        LU_2016mo(i:end,2) = LU_2016mo(i,2);
        
    end
    
    
    
else
    if strcmp(landusedata,'hansis')
      LU_2016 = csvread('Pongratz2016_GCP_meanPasture_peat.csv'); %Hansis GCP
    
    elseif strcmp(landusedata,'gcp')
        LU_2016 = csvread('GCPv1.3_historicalLU2016.csv');
        
    elseif strcmp(landusedata,'hough03')
        LU_2016 = csvread('LR_LU.csv'); %Rafelski 2009 high land use
        
    elseif strcmp(landusedata,'lr_low')
        LU_2016 = csvread('LR_LUex.csv'); %Rafelski 2009 low land use
    %elseif landusedata == 'const'
    
    elseif strcmp(landusedata,'CLM')
         CLMyear = [timeAuto; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(binomialDist), '.nc'),'time')];
    salinityAuto = [salinityAuto; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(binomialDist), '.nc'),'salinity')];

        
        
        
    end

    
    
luYear = LU_2016(1,1):(1/ts):LU_2016(end,1);
LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),luYear)).';
LU_2016mo(:,1) = luYear;
LU_2016mo(:,2) = LU_2016mo_0*d; % convert from PgC to ppm

end

    
        

        

% lu data, TgC/yr
%LU_2016 = csvread('GCPv1.3_historicalLU2016.csv');
%LU_2016 = csvread('Pongratz2016_GCP_meanPasture_peat.csv'); %Hansis GCP
%LU_2016 = csvread('LR_LU.csv'); %Rafelski 2009 high land use
%LU_2016 = csvread('LR_LUex.csv'); %Rafelski 2009 low land use
% luYear = LU_2016(1,1):(1/ts):LU_2016(end,1);
% LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),luYear)).';
% LU_2016mo(:,1) = luYear;
% LU_2016mo(:,2) = LU_2016mo_0*d; % convert from PgC to ppm




% to use a constant land use
% const = [LU_2016mo(1,1) LU_2016mo(1,2) ; LU_2016mo(end,1) LU_2016mo(end,2)];
% LU_2016mo(:,2) = (interp1(const(:,1),const(:,2),luYear)).';
% i = find(LU_2016mo(:,1) == 1920);
% LU_2016mo(i:end,2) = LU_2016mo(i,2);


% shorten datasets to match time frame of year vector
FF_start = find(FF_2016mo(:,1) == year(1));
FF_end = find(FF_2016mo(:,1) == year(end));
ff = FF_2016mo(FF_start:FF_end,:);
LU_start = find(LU_2016mo(:,1) == year(1));
LU_end = find(LU_2016mo(:,1) == year(end));
LU = LU_2016mo(LU_start:LU_end,:);

end