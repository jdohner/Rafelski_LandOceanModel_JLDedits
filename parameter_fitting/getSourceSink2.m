% getSourceSink2.m
%
% author: Lauren Rafelski, modified by Julia Dohner
% January 23, 2018
%
% reorganized version of getsourcesink_scale3 from LR code. Allows for
% forward model running to 2016
%
% instead of joining records at 1959, this version joins at latest
% timepoint

% LU in ppm/year
% ff in ppm/year


function [fas,ff,LU,LUex] = getSourceSink(year2);

load joos_hilda_2011.mat % ocean stuff    
    
% load data thru 2009, 2006, 2000 - ff is monthly, lu is annual

% both vectors begin in 1800
LU_2006 = csvread('landUse_1800-2006.csv');
LUex_2000 = csvread('landUseExtra_1800-2000.csv');
load fossilFuel_1751-2009.mat;
FF_2009 = ff1; % already monthly resolution

% shortening ff vector to begin at start_year (have data back thru 1700)
FF_start = find(FF_2009(:,1) == year2(1));
FF_2009 = FF_2009(FF_start:end,:);

% trimming start and end years
if year2(end) == 2016
    % shortening ff vector to begin at start_year
    FF_start = find(FF_2009(:,1) == year2(1));
    %FF_end = find(FF_2009(:,1) == year(end)); % buggy line - need data thru 2016
    %FF_2009 = FF_2009(FF_start:FF_end,:);
    FF_2009 = FF_2009(FF_start:end,:);
    ff = FF_2009;
else  
    FF_end = find(FF_2009(:,1) == year2(end)); % buggy line - need data thru 2016
    FF_2009 = FF_2009(1:FF_end,:);
    ff = FF_2009;
end

    
% interpolate to monthly
LU_2006mo(:,1) = (1800:(1/12):2006)';
LU_2006mo(:,2) = interp1(LU_2006(:,1),LU_2006(:,12),LU_2006mo(:,1));

LUex_2000mo(:,1) = (1800:(1/12):2000)';
LUex_2000mo(:,2) = interp1(LUex_2000(:,1),LUex_2000(:,2),LUex_2000mo(:,1));

% extend land use with last value thru to latest year
LU(:,1) = year2;
LU(1:length(LU_2006mo),2) = LU_2006mo(:,2);
LU(length(LU_2006mo)+1:end,2) = LU_2006mo(end,2);
    
% extend extratropical land use to end_year with zeros thru to latest year
LUex(:,1) = year2; 
LUex(1:length(LUex_2000mo),2) = LUex_2000mo(:,2);
LUex(length(LUex_2000mo)+1:end,2) = 0;


if year2(end) == 2016 % assuming either doing 2007 or 2016
    
    % load extended data thru 2016 - all annual
    FF_2016 = csvread('fossilFuel_1959-2016.csv'); % in gigatons/year
    % 1 ppm CO2 = 2.31 gton CO2
    d = 1/2.31; % gigaton to ppm conversion factor
    
    % this one isn't working for csvread, so reading in as text file
    %landUse_2016 = csvread('landUse_1959-2016.csv'); 
    fid = fopen('landUse_1959-2016.txt');
    C = textscan(fid,'%f %f', 'delimiter','\t');
    fclose(fid);
    LU_2016(:,1) = C{1};
    LU_2016(:,2) = C{2};
    
    % prep the 1959-2016 data vectors
    
    % fossil fuel 2016 to monthly
    month_2016 = 1959:(1/12):2016;
    FF_2016mo_0 = (interp1(FF_2016(:,1),FF_2016(:,2),month_2016)).';
    FF_2016mo(:,1) = month_2016;
    FF_2016mo(:,2) = FF_2016mo_0*d; %convert to ppm

    % land use 2016 to monthly
    LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),month_2016)).';
    LU_2016mo(:,1) = month_2016;
    LU_2016mo(:,2) = LU_2016mo_0*d;
    
    % extend (making full length vectors (1800-2016)) and patch

    clear ff; % get rid of dimension mismatch since looks like below pulls
    % from original data source (FF_2009) anyways
    % Houghton record begins in 1959, combine two records at 1959
    %i_2009 = find(FF_2009(:,1) == 2009); % same for ff, landuse, extraLU
    i_2009 = length(FF_2009);
    enddate_2009 = FF_2009(end,1);
    startdate_2016 = find(FF_2016mo(:,1) == enddate_2009);
    ff(:,1) = year2;
    ff(1:i_2009,2) = FF_2009(:,2);
    ff(i_2009+1:end,2) = FF_2016mo(startdate_2016+1:end,2);

    LU(:,1) = year2; 
    % just extend LU from last value, don't use new data
    LU(1:length(LU_2006mo),2) = LU_2006mo(:,2);
    LU(length(LU_2006mo)+1:end,2) = LU_2006mo(end,2);
    
    
%     LU(1:i_2009-1,2) = LU_2006mo(1:i_2009-1,2);
%     LU(i_2009:end,2) = LU_2016mo(:,2);

end

