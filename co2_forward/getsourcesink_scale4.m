% Get relevant land use and Joos model results
%
% Updates
%
% 4/13/07 - Lauren Elmegreen - use land use change 
% emissions that are extrapolated to the present by keeping values 
% constant at the last value, instead of increasing by 1.4%
% 10/25/12 - Lauren Rafelski - added annotations
% Oct 24, 2017 - JLD removed Aoc from output
%
% This file is for the record between 1800 and 2006+10/12

function [ff,LU,LUex] = getsourcesink_scale4(predict,year);

if predict == 1

    % load data thru 2009, 2006, 2000 - ff is monthly, lu is annual

    % both vectors begin in 1800
    LU_2006 = csvread('landUse_1800-2006.csv');
    LUex_2000 = csvread('landUseExtra_1800-2000.csv');
    load fossilFuel_1751-2009.mat;
    FF_2009 = ff1; % already monthly resolution
    % shortening ff vector to begin at start_year
    FF_start = find(FF_2009(:,1) == year(1));
    FF_2009 = FF_2009(FF_start:end,:);

    % load extended data thru 2016 - all annual
    FF_2016 = csvread('fossilFuel_1959-2016.csv');
    % this one isn't working for csvread, so reading in as text file
    %landUse_2016 = csvread('landUse_1959-2016.csv'); 
    fid = fopen('landUse_1959-2016.txt');
    C = textscan(fid,'%f %f', 'delimiter','\t');
    fclose(fid);
    LU_2016(:,1) = C{1};
    LU_2016(:,2) = C{2};


    % interpolate to monthly

    % landuse 2006 to monthly
    month_2006 = 1800:(1/12):2006;
    LU_2006mo = interp1(LU_2006(:,1),LU_2006(:,12),month_2006);
    LU_2006mo(:,1) = month_2006;
    LU_2006mo(:,2) = LU_2006mo; % value in ppm

    % extratrop landuse 2000 to monthly
    month_2000 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
    LUex_2000mo = interp1(LUex_2000(:,1),LUex_2000(:,2),month_2000);
    LUex_2000mo(:,1) = month_2000;
    LUex_2000mo(:,2) = LUex_2000mo; % value in ppm

    % fossil fuel 2009 already in monthly

    % fossil fuel 2016 to monthly
    month_2016 = 1959:(1/12):2016;
    FF_2016mo = (interp1(FF_2016(:,1),FF_2016(:,2),month_2016)).';
    FF_2016mo(:,1) = month_2016;
    FF_2016mo(:,2) = FF_2016mo;

    % land use 2016 to monthly
    LU_2016mo = (interp1(LU_2016(:,1),LU_2016(:,2),month_2016)).';
    LU_2016mo(:,1) = month_2016;
    LU_2016mo(:,2) = LU_2016mo;



    % extend (making full length vectors (1800-2016)) and patch

    % Houghton record begins in 1959, combine two records at 1959
    i_1959 = find(FF_2009(:,1) == 1959); % same for ff, landuse, extraLU
    ff(:,1) = year;
    ff(1:i_1959-1,2) = FF_2009(1:i_1959-1,2);
    ff(i_1959:end,2) = FF_2016mo(:,2);

    LU(:,1) = year; 
    LU(1:i_1959-1,2) = LU_2006mo(1:i_1959-1,2);
    LU(i_1959:end,2) = LU_2016mo(:,2);

    % don't have updated extratropical land use data, so just extend from 2000
    % onward using 0
    LUex(:,1) = year; 
    LUex(1:length(LUex_2000mo),2) = LUex_2000mo(:,2);
    LUex(length(LUex_2000mo)+1:end,2) = 0;



else % if predict == 0 (diagnostic case, not forward projecting), LR code
    
    load joos_hilda_2011.mat
    landnowppm_JLD = csvread('landnowppm_JLD.csv');
    extratrop_landppm_JLD = csvread('extratrop_landppm_JLD.csv');


    % interpolate land use changes to monthly resolution

    month_2006 = 1800:(1/12):2006;
    landmonth_2006 = interp1(landnowppm_JLD(:,1),landnowppm_JLD(:,12),month_2006);

    month_2000 = 1800:(1/12):2000; % extratrop record ends in 2000, will extend with 0's same way LR did
    extralandmonth_2000 = interp1(extratrop_landppm_JLD(:,1),extratrop_landppm_JLD(:,2),month_2000);

    LU_2006mo(:,1) = month_2006;
    LU_2006mo(:,2) = landmonth_2006; % value in ppm

    LUex_2000mo(:,1) = month_2000;
    LUex_2000mo(:,2) = extralandmonth_2000; % value in ppm

    last_ind_2006 = length(month_2006);
    last_ind_2000 = length(month_2000);

    %Extend extratropical emissions by assuming emissions are zero
    %extratrop_landmo(:,1) = landusemo(:,1); % extend time column
    LUex_2000mo(last_ind_2000+1:last_ind_2006,1) = LU_2006mo(last_ind_2000+1:last_ind_2006,1);
    LUex_2000mo(last_ind_2000+1:last_ind_2006,2) = 0;

    LU = LU_2006mo;
    LUex = LUex_2000mo;
    ff = ff1;
    % other variables should just be loaded and passed

end

%%%%
% getsourcesink outputs end at 2006
% Extend land use record by making recent emissions equal to last
% record
j = length(LU);
%k = length(year);
%extendYears = year(1,j+1:end);
%landusemo(:,1) = [landusemo(:,1), extendYears];

%LU(j+1:k,1) = year(1,j+1:k);
LU(j+1:end,2) = LU(j,2);

%Extend extratropical emissions by assuming emissions are zero
i = length(LUex);
%LUex(i+1:k,1) = LU(i+1:k,1);
LUex(i+1:end,2) = 0;


% TODO: the following (extending records) can all happen in 
% getsourcesink_scale3, but wait to
% change until can get the code running and working