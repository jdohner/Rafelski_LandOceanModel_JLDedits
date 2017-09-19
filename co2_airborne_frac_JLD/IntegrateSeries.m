% file integrate_series_trap2.m
%
% author Lauren Rafelski, modified by Julia Dohner
% 
% brief integrate_series_trap2 is a function that integrates data based on
% the trapezoid rule. Time points and data points should be in the same
% input array, with time points in the column "timecol" and data points in
% the column "numcol". dt is the number of points per year.
%
% TODO: rename as IntegrateSeries.m

function [func_int] = IntegrateSeries(func,timecol,numcol,dt)

func_int = zeros(length(func),numcol);
func_int(:,1) = func(:,timecol);

for i = 2:length(func(:,timecol))
    func_int(i,2) = func_int(i-1,2) + 0.5*(func(i-1,numcol)+func(i,numcol))/dt;
end