% tempRecord.mat
%
% author: Lauren Rafelski, modified by Julia Dohner
% Jan 24, 2018
%
% puts all the temperature record code into separate function to clean up
% main driver
%
% CRUTEM4-gl data from https://crudata.uea.ac.uk/cru/data/temperature/
% Climatic Research Unit
%
% notes: could clean up to make it run through full history, then cut at
% the end based on end_year
% also worth noting that CRU data has much better resolution- should use
% throughout or just append to after 2009?

function [temp_anom, T0] = tempRecord2(start_year,end_year,dt)

load landwt_T_2011.mat

CRU_data = csvread('CRUTEM4-gl.csv');
    
% processing of CRUTEM4 file - comes in weird format
% taking just rows with temp anomalies and cut off year values
CRU_startYr = CRU_data(1,1); % starts 1850
CRU_endYr = CRU_data(end,1) + (11/12); % ends Dec 2017
CRU_year = (CRU_startYr:(1/12):CRU_endYr)';
CRU_temp = CRU_data(1:2:end,2:13); % every other row for cols 2 to 13
CRU_temp = reshape(CRU_temp',[],1); % reshape to column vector

% temp_anom is 1850-Dec 2017 (end CRU data)
temp_full(:,1) = 1850:dt:CRU_endYr;
% add in landtglob, which begins at 1850.5
i = find(temp_full(:,1) == landtglob(1,1));
temp_full(1:i-1,2) = landtglob(1,2); % set first values
% add landtglob up until 2009+7/12
j = find(floor(100*temp_full(:,1)) == floor(100*(2009+(7/12))));
k = find(floor(100*landtglob(:,1)) == floor(100*(2009+(7/12)))); 
temp_full(i:j,2) = landtglob(1:k,2);  %%
l = find(floor(100*CRU_year) == floor(100*landtglob(1910,1))); % find 2009+7/12
%m = find(CRU_year == end);
temp_full(1917:end,2) = CRU_temp(l+1:end);

% temp_anom record spans 1850-dec 2017
% cut short according to start and end years
i1 = find(floor(100*temp_full(:,1)) == floor(100*start_year));
j1 = find(floor(100*temp_full(:,1)) == floor(100*end_year));
temp_anom = temp_full(i1:j1,:);

T0 = temp_full(1,2); % for some reason want first point at 1850