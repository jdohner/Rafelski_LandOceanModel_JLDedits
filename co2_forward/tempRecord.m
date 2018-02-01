% tempRecord.mat
%
% author: Lauren Rafelski, modified by Julia Dohner
% Jan 24, 2018
%
% puts all the temperature record code into separate function to clean up
% main driver

function [temp_anom, T0] = tempRecord(tRecord,globalRecord,dt, end_year);

if end_year == 2016
    globT_2017 = csvread('HadCRUT4_2017.csv');
end


% do a moving boxcar average of the land temperature: 1 year average
%l_boxcar(func,boxlength,dt,starttime,endtime,datecol,numcol)
[avg_temp] = l_boxcar(tRecord,1,12,1,2483,1,2); 
% at the moment, length(year_ocean) is too long for the tland4 record
%[avg_temp] = l_boxcar(tempRecord,1,12,1,length(year_ocean),1,2); 

% THIS IS REALLY SUSPECT - looks like globalRecord is just repeated
k = find(avg_temp(:,1) == globalRecord(1,1));
temp_anom(1:k-1,1) = avg_temp(1:k-1,1); % fill time column of temp_anom
% fill with first value from globalRecord
temp_anom(1:k-1,2) = globalRecord(1,2); %globalRecord are temp anomalies
j = find(floor(100*globalRecord(:,1)) == floor(100*end_year));
temp_anom(k:j+k-1,1) = globalRecord(1:j,1); % globalRecord begins 1850.5
temp_anom(k:j+k-1,2) = globalRecord(1:j,2);

clear avg_temp;

T0 = temp_anom(1,2);