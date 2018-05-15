% New subroutine friendly version of joosHILDA - variable temp
% 1/29/08: changed temperature to new record


% addpath(genpath(...
%     '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data'));
% addpath(genpath(...
%     '/Users/juliadohner/Documents/MATLAB/JLDedits_Rafelski_LandOceanModel/co2_forward/co2_forward_data_2016'));

function [fas] = joosHILDAvarT_scale_JLD(start_year,end_year,ts,ff);

ts = 12;
start_yearOcean = 1800;
%end_year = 2010;%2005.5;
Aoc_HvarT = 3.62E14;
c = 1.722E17;
h = 75;
kg = 1/9.06;
meanSST = 18.35;
yearOcean = start_yearOcean:(1/ts):end_year;

[dtdelpCO2a,dpCO2a,year,dt,~] = getObservedCO2_2(ts,start_yearOcean,end_year);

load newsst_2011.mat;
hadcrut = csvread('HadCRUT4.csv');
% interpolate to monthly
sstmo(:,1) = yearOcean;
sstmo(:,2) = interp1(hadcrut(1:2:end,1),hadcrut(1:2:end,14),yearOcean);
i = find(sstmo(:,1) == 1850);
sstmo(1:i-1,2) = 0; % fill in values before 1850
sstmo(:,2) = sstmo(:,2) + meanSST;
T = sstmo;

% don't need to smooth
% [T] = l_boxcar(sstmo,1,12,1,length(sstmo),1,2);
% T(1:6,1) = [1800:1/12:(1800.5-(1/12))];
% T(1:6,2) = 17.967;

% Response function
for n = 1:length(year)
     t(n,1) = year(n) - year(1);
     r(n,1) = t(n,1);
    %Calculate response function based on HILDA equation
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

[fas_HvarT] = joos_varT2_water(yearOcean,dpCO2a,c,h,kg,T,Aoc_HvarT,r,dt);

% Calculate land flux using fossil fuel sources and ocean sink (ocean
% flux*area of ocean)


for p = 1:(length(year)-1)
    q = find(ff1(:,1) == year(1,p));
    landflux_HvarT(p,1) = year(p);
    landflux_HvarT(p,2) = dtdelpCO2a(p+((1800-1640)*ts),2) - ff1(q,2) + fas_HvarT(p,2)*Aoc_HvarT; 
end

% shorten to start and end years

%x = 0*year;

% figure
% plot(ff1(:,1),ff1(:,2),'-k',dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r',fas_HvarT(:,1),-Aoc_HvarT*fas_HvarT(:,2),'-b',landflux_HvarT(:,1),landflux_HvarT(:,2),'-g',year(1,:),x,'--k')
% axis([1800 2010 -10 10])
% legend('fossil fuel','atmosphere','ocean','land','Location','SouthWest')
% title('Sources and sinks from Joos HILDA model with variable T ')
% xlabel('Year ')
% ylabel('ppm/year  Positive = source, negative = sink ')