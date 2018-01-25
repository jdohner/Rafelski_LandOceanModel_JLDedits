% 5/9/07 - change so that it always puts the date in column 1 and the data in column 2

function [avg_func] = l_boxcar(func,boxlength,dt,starttime,endtime,datecol,numcol)


for i = (starttime+(boxlength/2)*dt):(endtime-(boxlength/2)*dt)
    avg_func(i,1) = func(i,datecol);
    avg_func(i,2) = sum(func(i-(boxlength/2)*dt:i+(boxlength/2)*dt,numcol))/(boxlength*dt+1);
end

% this boxcar function is buggy and i doesn't start at 1, so first few
% values are left unfilled

% to fill the first values:
idx = find(avg_func~=0, 1, 'first'); % find index of first non-zero element
avg_func(1:idx-1,2) = avg_func(idx,2);