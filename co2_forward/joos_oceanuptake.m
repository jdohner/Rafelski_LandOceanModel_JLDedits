% joos_oceanuptake
%
% parallel to joos_general_fast_annotate2.n
% calculates fas for all timepoints through to present (2017)
% author: lr adapted by jld
% april 4, 2018

% Pulse response model using response function "r"

%[fas,dpCO2s] = joos_general_fast_annotate2(year,dpCO2a,c,h,kg,T,Aoc,r,dt);
function [fas,dpCO2s]= joos_oceanuptake(year,c,h,kg,T,Aoc,dt)
%load joos_hilda_2011.mat % ocean stuff
end_year = year(end,1);
ts = 12;
[~,dpCO2a,~] = MLOinterpolate_increment2_recent(ts,1800,end_year); 


% response function to calculate ocean uptake
% Response function to calculate ocean uptake
for n = 1:length(year)
     t(n,1) = year(n) - year(1);
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
    
% rest of ocean uptake calculation stuff 
    
dpCO2s = zeros(length(dpCO2a),2); % dissolved CO2
dpCO2s(:,1) = dpCO2a(:,1);
fas = zeros(length(year),2);
integral = zeros(length(year),length(year));
delDIC = zeros(length(year),2);

for m = 1:(length(year)-1)
    %Calculate flux 
   fas(m,1) = year(m);
   fas(m,2) = (kg/Aoc)*(dpCO2a(m,2) - dpCO2s(m,2)); % air-sea flux of CO2
   
    w = conv(fas(1:m,2),r(1:m,2)); % convolve the air-sea flux and the pulse response function, as in Joos 1996

    % Calculate delDIC 
    delDIC(m+1,1) = year(m+1);
   
    delDIC(m+1,2) = (c/h)*w(m)*dt; % change in DIC

    %Calculate dpCO2s from DIC - from Joos 1996
    dpCO2s(m+1,2) = (1.5568 - (1.3993E-2)*T)*delDIC(m+1,2) + (7.4706-0.20207*T)*10^(-3)*...
        (delDIC(m+1,2))^2 - (1.2748-0.12015*T)*10^(-5)*(delDIC(m+1,2))^3 + (2.4491-0.12639*T)...
        *10^(-7)*(delDIC(m+1,2))^4 - (1.5468-0.15326*T)*10^(-10)*(delDIC(m+1,2))^5;
    
end

% calculate the flux for the last time point
fas(length(year),1) = year(length(year));
fas(length(year),2) = (kg/Aoc)*(dpCO2a(length(year),2) - dpCO2s(length(year),2));

zero1 = find(delDIC(:,1) == 0);
delDIC(zero1,2) = NaN;
delDIC(zero1,1) = NaN;

zero2 = find(dpCO2s(:,2) == 0);
dpCO2s(zero2,2)  = NaN;
