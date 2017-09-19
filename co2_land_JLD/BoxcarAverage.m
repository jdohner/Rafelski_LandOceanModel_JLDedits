% file BoxcarAverage.m
% formerly "l_boxcar.m"
%
% author Lauren Rafelski, modified by Julia Dohner
%
% brief The function BoxcarAverage gives a moving boxcar average of a 
% temperature record.
%
% "Boxcar averaging is a signal smoothing technique that assumes the
% average of a small numer of adjacent points to be a better measure of
% signal than any of the individual points." (Source: 
% https://sites.google.com/a/cord.edu/labview-for-analytical-chemistry/ ...
% home/simulations/boxcar-averaging)
%
% inputs tland4 (temperature record), boxlength (in units of years), dt 
% (number of points per year), starttime, endtime, datecol (to indicade
% which column the date should appear in), and numcol (to indicate which 
% column the data should appear in).
% 
% outputs avg_temp (first column of avg_temp gives the date, second column
% gives the moving average of the land temperature).
% 
% changes by LR:
% 5/9/07 - change so that it always puts the date in column 1 and the data 
% in column 2
%
% TODO: can edit this so that always put date in column 1 and data in
% column 2 without having to take dateCol and numCol as inputs

function [avg_temp] = BoxcarAverage(landTempRecord,boxLength,pointsPerYear,startTime,endTime,dateCol,numCol)

window = boxLength*pointsPerYear; % window of data within which taking average
numRows = (endTime-(window/2)) - (startTime+(window/2));
avg_temp = NaN(numRows, 2); % allocate space for avg_temp matrix before loop

% taking average of point exactly in middle of window (hence +/- window/2)
for i = (startTime+(window/2)):(endTime-(window/2))
    avg_temp(i,1) = landTempRecord(i,dateCol); % add date to avg_temp array
    avg_temp(i,2) = sum(landTempRecord(i-(window/2):i+(window/2),numCol))/(window+1);
end
