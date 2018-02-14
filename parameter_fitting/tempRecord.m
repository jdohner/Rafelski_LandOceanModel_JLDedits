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

function [temp_anom, T0] = tempRecord(temp_early,temp_recent,dt, year, end_year);

if end_year == 2016
    clear temp_recent;
    CRU_data = csvread('CRUTEM4-gl.csv');

    % processing of CRUTEM4 file
    % taking just rows with temp anomalies and cut off year values
    CRU_startYr = CRU_data(1,1);
    CRU_endYr = CRU_data(end,1) + (11/12); % ends Dec 2017
    CRU_year = (CRU_startYr:(1/12):CRU_endYr)';
    CRU_temp = CRU_data(1:2:end,2:13); % every other row for cols 2 to 13
    CRU_temp = reshape(CRU_temp',[],1); % reshape to column vector
    temp_recent(:,1) = CRU_year;
    temp_recent(:,2) = CRU_temp;

end

% do a moving boxcar average of the land temperature: 1 year average
%l_boxcar(func,boxlength,dt,starttime,endtime,datecol,numcol)
[avg_temp] = l_boxcar(temp_early,1,12,1,2483,1,2); 
% at the moment, length(year_ocean) is too long for the tland4 record
%[avg_temp] = l_boxcar(tempRecord,1,12,1,length(year_ocean),1,2); 

% THIS IS REALLY SUSPECT - looks like globalRecord is just repeated
k = find(avg_temp(:,1) == temp_recent(1,1));
temp_anom(1:k-1,1) = avg_temp(1:k-1,1); % fill time column of temp_anom
% fill with first value from globalRecord
temp_anom(1:k-1,2) = temp_recent(1,2); %globalRecord are temp anomalies
j = find(floor(100*temp_recent(:,1)) == floor(100*end_year));
temp_anom(k:j+k-1,1) = temp_recent(1:j,1); % globalRecord begins 1850.5
temp_anom(k:j+k-1,2) = temp_recent(1:j,2);

clear avg_temp;

T0 = temp_anom(1,2);