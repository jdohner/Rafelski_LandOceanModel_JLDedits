% file OceanPulseResponse.m
% formerly "joos_general_fast_annotate2.m"
%
% author Lauren Rafelski, modified by Julia Dohner
% 
% inputs year,dpCO2a,c,h,kg,T,Aoc,r,dt
% outputs [airSeaFlux,dpCO2s]
% 
% brief The function joos_general_fast_annotate2 is the ocean uptake 
% pulse response model "r" from Joos (1996) defined in section A.2.2. 
% (page 415) of the paper.
% 
% changes by LR
% 7/15/09: add "dpCO2s" to output

function [airSeaFlux,dpCO2s]= OceanPulseResponse(year,dpCO2a,c,h,kg,T,Aoc,r,dt)
%input dpCO2a is change in atmospheric pCO2

%initialize dpCO2s, which is change in dissolved pCO2
%set dpCO2s to be same dimensions as dpCO2a
dpCO2s = zeros(length(dpCO2a),2); 
% fill first column of dpCO2s with same values as dpCO2a
dpCO2s(:,1) = dpCO2a(:,1); 

airSeaFlux = zeros(length(year),2);
airSeaFlux(:,1) = year;
integral = zeros(length(year),length(year));
delDIC = zeros(length(year),2);

% Calculate flux
for m = 1:(length(year)-1)
    
   % airSeaFlux(m,1) = year(m); JLD
   % Seems like there's an issue here: about here, dpCO2a and dpCO2s are the
   % same. I don't see how they'd be any different given the above code.
   % Update: after some debugging, it looks like dpCO2a has values in the
   % second column while dpCO2s only has zeros in the second column. Is
   % this scientifically correct? The ocean's change in co2 is 0?
   airSeaFlux(m,2) = (kg/Aoc)*(dpCO2a(m,2) - dpCO2s(m,2)); % air-sea flux of CO2
   
   % convolve the air-sea flux and the pulse response function, as in Joos 1996
   w = conv(airSeaFlux(1:m,2),r(1:m,2)); 
    

   % Calculate delDIC 
   delDIC(m+1,1) = year(m+1);
   
   % Perhaps equation (3) from Joos (1996) page 401
   delDIC(m+1,2) = (c/h)*w(m)*dt; % change in DIC

   % Calculate dpCO2s from DIC
   % Equation (6b) in Joos (1996) page 402
   dpCO2s(m+1,2) = (1.5568 - (1.3993E-2)*T)*delDIC(m+1,2) + (7.4706-0.20207*T)*10^(-3)*...
       (delDIC(m+1,2))^2 - (1.2748-0.12015*T)*10^(-5)*(delDIC(m+1,2))^3 + (2.4491-0.12639*T)...
       *10^(-7)*(delDIC(m+1,2))^4 - (1.5468-0.15326*T)*10^(-10)*(delDIC(m+1,2))^5;
    
end

% calculate the flux for the last time point
airSeaFlux(end,1) = year(end);
airSeaFlux(end,2) = (kg/Aoc)*(dpCO2a(end,2) - dpCO2s(end,2));

% find index of zero in delDIC matrix
zero1 = find(delDIC(:,1) == 0);
% change elements in row where value 0 appears to be NaN's
delDIC(zero1,2) = NaN;
delDIC(zero1,1) = NaN;

% find index of zero in dpCO2s matrix
zero2 = dpCO2s(:,2) == 0;
% change element at index zero2 to be NaN
dpCO2s(zero2,2)  = NaN;