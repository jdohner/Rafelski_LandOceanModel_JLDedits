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


% do a moving boxcar average of the land temperature: 1 year average
avg_temp = l_boxcar(temp_early,1,12,1,2483,1,2); 
avg_temp(1:6,2) = avg_temp(7,2); % make the first 6 points 

% THIS IS REALLY SUSPECT - looks like globalRecord is just repeated
k = find(avg_temp(:,1) == temp_recent(1,1));
temp_anom(1:k-1,1) = avg_temp(1:k-1,1); % fill time column of temp_anom
% fill with first value from globalRecord
temp_anom(1:k-1,2) = temp_recent(1,2); %globalRecord are temp anomalies
j = find(floor(100*temp_recent(:,1)) == floor(100*end_year));
temp_anom(k:j+k-1,1) = temp_recent(1:j,1); % globalRecord begins 1850.5
temp_anom(k:j+k-1,2) = temp_recent(1:j,2);

if end_year >= 2009+(7/12)
    %clear temp_recent;
    CRU_data = csvread('CRUTEM4-gl.csv');

    % processing of CRUTEM4 file - comes in weird format
    % taking just rows with temp anomalies and cut off year values
    CRU_startYr = CRU_data(1,1); % starts 1850
    CRU_endYr = CRU_data(end,1) + (11/12); % ends Dec 2017
    CRU_year = (CRU_startYr:(1/12):CRU_endYr)';
    CRU_temp = CRU_data(1:2:end,2:13); % every other row for cols 2 to 13
    CRU_temp = reshape(CRU_temp',[],1); % reshape to column vector
    
    i = find(floor(100*CRU_year) == floor(100*temp_early(end,1))); % deals with last decimal digit in repeating
    j = length(CRU_year)-i; % how many additional values adding to temp_recent
    temp_extend(:,1) = [temp_anom(:,1) ; CRU_year(i:end)];
    
    temp_recent(i:j,1) = CRU_year(i:end);
    temp_recent(:,2) = CRU_temp(i:end);
    % CRU data becomes the data for 1850-2017.9

end

clear avg_temp;

T0 = temp_anom(1,2);