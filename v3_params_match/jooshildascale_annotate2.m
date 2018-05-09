%% January 10, 2011: extended data to 2010



function [fas,sstAnom] = jooshildascale_annotate2(start_year,end_year,ts,ff,varSST,Tconst);

ts = 12; % number of data points/year
start_yearOcean = 1800;
%end_year = 2009+(7/12);
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996



%[dtdelpCO2a,dpCO2a,year,dt] = MLOinterpolate_increment2(ts,start_year,end_year); % get atmospheric CO2 record
[~,dpCO2a,yearOcean,dt,~] = getObservedCO2_2(ts,start_yearOcean,end_year);

%[ff1] = load_fossil2(ts); % get fossil fuel emissions


% Response function to calculate ocean uptake
for n = 1:length(yearOcean)
     t(n,1) = yearOcean(n) - yearOcean(1);
     r(n,1) = t(n,1);
    %Calculate response function based on HILDA equation in Joos 1996
     if t(n,1) == 0
         r(n,2) = 1;
     elseif t(n,1) <= 2
         r(n,2)= (1/0.95873)*(0.12935+0.21898*exp(-t(n,1)/0.034569)+0.17003*exp(-t(n,1)/0.26936)...
             +0.24071*exp(-t(n,1)/0.96083)+0.24093*exp(-t(n,1)/4.9792));
     else
         r(n,2) = (1/0.95873)*(0.022936+0.24278*exp(-t(n,1)/1.2679)+0.13963*exp(-t(n,1)/5.2528)...
                +0.089318*exp(-t(n,1)/18.601)+0.037820*exp(-t(n,1)/68.736)...
                +0.035549*exp(-t(n,1)/232.3));
     end
end

if varSST == 0
    T = Tconst;
end




%% calculate ocean uptake
[fas,dpCO2s,sstAnom] = joos_general_fast_annotate2(yearOcean,dpCO2a,c,h,kg,Tconst,Aoc,r,dt,varSST); 

i3 = find(fas(:,1) >= start_year,1);
j3 = find(fas(:,1) >= end_year,1);
fas = fas(i3:j3,:);



end
